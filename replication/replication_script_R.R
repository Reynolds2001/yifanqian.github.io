############################################################
# Replication script: Russia opposition network appendix

install_missing <- TRUE
required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "purrr", "stringr",
  "tibble", "sandwich"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  if (install_missing) {
    install.packages(missing_packages)
  } else {
    stop("Please install these packages before running the script: ",
         paste(missing_packages, collapse = ", "))
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(sandwich)
})

############################
# 1. File paths
############################

nodes_path <- "network_nodes.csv"
edges_path <- "network_edges.csv"
output_dir <- "replication_output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

############################
# 2. Helper functions
############################

read_network_table <- function(path, sheet = 1) {
  # Some uploaded files have a .csv extension but are actually Excel workbooks.
  # Try Excel first, then CSV.
  x <- tryCatch(
    suppressMessages(readxl::read_excel(path, sheet = sheet)),
    error = function(e) NULL
  )

  if (is.null(x)) {
    x <- tryCatch(
      suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE)),
      error = function(e) NULL
    )
  }

  if (is.null(x)) {
    stop("Could not read file: ", path,
         "\nIf this is an Excel workbook with a .csv extension, rename it to .xlsx or use the original .xlsx file.")
  }

  as_tibble(x)
}

pop_sd <- function(x) {
  # Population standard deviation, not R's default sample SD.
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  sqrt(mean((x - mean(x))^2))
}

make_dyads_one_year <- function(.x, .y) {
  actors <- sort(unique(.x$actor))
  cmb <- t(utils::combn(actors, 2))
  tibble(
    year = as.integer(.y$year[[1]]),
    actor_i = cmb[, 1],
    actor_j = cmb[, 2]
  )
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 4)))
}

stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

fit_extract <- function(data, formula, term, se_type = c("ols", "hc1", "cluster"), cluster_var = NULL) {
  se_type <- match.arg(se_type)
  dat <- as.data.frame(data)
  dat$.row_id__ <- seq_len(nrow(dat))

  model <- stats::lm(formula, data = dat)
  coef_vec <- stats::coef(model)

  if (!term %in% names(coef_vec)) {
    stop("Term not found in model coefficients: ", term,
         "\nAvailable terms: ", paste(names(coef_vec), collapse = ", "))
  }

  if (se_type == "ols") {
    vc <- stats::vcov(model)
    df <- stats::df.residual(model)
    inference_rule <- "Conventional OLS SE"
  }

  if (se_type == "hc1") {
    vc <- sandwich::vcovHC(model, type = "HC1")
    df <- stats::df.residual(model)
    inference_rule <- "HC1"
  }

  if (se_type == "cluster") {
    if (is.null(cluster_var)) stop("cluster_var must be provided when se_type = 'cluster'.")
    mf <- stats::model.frame(model)
    model_rows <- as.integer(rownames(mf))
    cluster <- dat[[cluster_var]][model_rows]
    vc <- sandwich::vcovCL(model, cluster = cluster, type = "HC1")
    G <- length(unique(cluster[!is.na(cluster)]))
    df <- G - 1
    inference_rule <- dplyr::case_when(
      cluster_var == "actor" ~ "Actor-clustered",
      cluster_var == "dyad_id" ~ "Dyad-clustered",
      TRUE ~ paste0("Clustered: ", cluster_var)
    )
  }

  beta <- unname(coef_vec[term])
  se <- sqrt(diag(vc))[term]
  t_stat <- beta / se
  p_value <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)

  tibble(
    term = term,
    coefficient = beta,
    se = unname(se),
    t_stat = unname(t_stat),
    p_value = unname(p_value),
    df = df,
    n = stats::nobs(model),
    inference_rule = inference_rule
  )
}

add_display_cols <- function(x) {
  x %>%
    mutate(
      coef_display = paste0(fmt_num(coefficient), stars(p_value)),
      se_display = if_else(is.na(se), "—", paste0("(", fmt_num(se), ")")),
      p_display = fmt_p(p_value)
    )
}

run_spec <- function(table, model, outcome, specification, predictor_label,
                     data, formula, term, se_type, cluster_var = NULL) {
  fit_extract(data, formula, term, se_type, cluster_var) %>%
    mutate(
      table = table,
      model = model,
      outcome = outcome,
      specification = specification,
      predictor_label = predictor_label,
      year_fe = "Yes"
    ) %>%
    add_display_cols() %>%
    select(table, model, outcome, specification, predictor_label,
           coefficient, se, p_value, t_stat, df, n, inference_rule,
           coef_display, se_display, p_display, year_fe)
}

# Permutation inference for H3/H4.
# Logic: within each election cycle, randomly reassign actor-level adverse-treatment
# positions across actors, rebuild the dyad-level both-high or both-low indicator,
# re-estimate the same year-FE LPM, and compare the observed coefficient to the
# randomization distribution.
permute_adverse_test <- function(table, model, outcome, specification, predictor_label,
                                 data, actor_data, observed_var,
                                 mode = c("high", "low"),
                                 B = 9999, seed = 42) {
  mode <- match.arg(mode)

  obs_formula <- stats::as.formula(paste0(outcome, " ~ ", observed_var, " + factor(year)"))
  obs_fit <- stats::lm(obs_formula, data = data)
  obs_beta <- unname(stats::coef(obs_fit)[observed_var])

  actor_base <- actor_data %>%
    select(year, actor, adverse_treatment)

  dyad_base <- data %>%
    select(year, actor_i, actor_j, all_of(outcome))

  set.seed(seed)

  perm_betas <- replicate(B, {
    perm_actor <- actor_base %>%
      group_by(year) %>%
      mutate(adverse_perm = sample(adverse_treatment, size = n(), replace = FALSE)) %>%
      ungroup() %>%
      select(year, actor, adverse_perm)

    perm_dyads <- dyad_base %>%
      left_join(
        perm_actor %>% rename(actor_i = actor, adverse_i = adverse_perm),
        by = c("year", "actor_i")
      ) %>%
      left_join(
        perm_actor %>% rename(actor_j = actor, adverse_j = adverse_perm),
        by = c("year", "actor_j")
      )

    if (mode == "high") {
      perm_dyads <- perm_dyads %>%
        mutate(perm_predictor = as.integer(adverse_i > 0 & adverse_j > 0))
    } else {
      perm_dyads <- perm_dyads %>%
        mutate(perm_predictor = as.integer(adverse_i <= 0 & adverse_j <= 0))
    }

    perm_fit <- stats::lm(
      stats::as.formula(paste0(outcome, " ~ perm_predictor + factor(year)")),
      data = perm_dyads
    )

    unname(stats::coef(perm_fit)["perm_predictor"])
  })

  # Two-sided randomization p-value, with +1 correction so p is never exactly zero.
  perm_p <- (1 + sum(abs(perm_betas) >= abs(obs_beta), na.rm = TRUE)) / (B + 1)

  tibble(
    table = table,
    model = model,
    outcome = stringr::str_to_sentence(stringr::str_replace_all(outcome, "_", " ")),
    specification = specification,
    predictor_label = predictor_label,
    coefficient = obs_beta,
    se = NA_real_,
    p_value = perm_p,
    t_stat = NA_real_,
    df = NA_real_,
    n = nrow(data),
    inference_rule = "Permutation",
    year_fe = "Yes"
  ) %>%
    add_display_cols() %>%
    select(table, model, outcome, specification, predictor_label,
           coefficient, se, p_value, t_stat, df, n, inference_rule,
           coef_display, se_display, p_display, year_fe)
}

write_table <- function(x, filename) {
  readr::write_csv(x, file.path(output_dir, filename))
}

############################
# 3. Read and clean input data
############################

nodes_raw <- read_network_table(nodes_path, sheet = 1)
edges_raw <- read_network_table(edges_path, sheet = 1)

nodes <- nodes_raw %>%
  filter(!is.na(year), !is.na(actor)) %>%
  mutate(
    year = as.integer(year),
    actor = as.character(actor),
    system_status = as.character(system_status),
    parliamentary = as.integer(parliamentary),
    opposition_actor = as.integer(opposition_actor)
  )

edges <- edges_raw %>%
  filter(!is.na(year), !is.na(source), !is.na(target)) %>%
  mutate(
    year = as.integer(year),
    source = as.character(source),
    target = as.character(target),
    edge_type = as.character(edge_type),
    sign = as.integer(sign),
    weight = as.numeric(weight),
    directed = as.logical(directed),
    evidence = as.character(evidence),
    source_url = as.character(source_url)
  )

cycle_order <- sort(unique(nodes$year))
prev_lookup <- tibble(
  year = cycle_order[-1],
  prev_year = cycle_order[-length(cycle_order)]
)

node_status <- nodes %>%
  select(year, actor, system_status, opposition_actor)

opposition_nodes <- nodes %>%
  filter(opposition_actor == 1)

############################
# 4. Dataset scope and edge summaries
############################

scope_by_cycle <- opposition_nodes %>%
  group_by(year) %>%
  summarise(
    opposition_actors = n_distinct(actor),
    opposition_dyads = opposition_actors * (opposition_actors - 1) / 2,
    .groups = "drop"
  )

edge_type_counts <- edges %>%
  count(edge_type, sign, name = "n") %>%
  arrange(sign, edge_type)

write_table(scope_by_cycle, "Table_A3_dataset_scope_by_cycle.csv")
write_table(edge_type_counts, "edge_type_counts.csv")

############################
# 5. Regime access and pressure construction
############################

regime_treatment <- edges %>%
  filter(source == "United Russia") %>%
  group_by(year, actor = target) %>%
  summarise(
    access = sum(if_else(sign == 1, weight, 0), na.rm = TRUE),
    pressure = sum(if_else(sign == -1, weight, 0), na.rm = TRUE),
    .groups = "drop"
  )

actor_data <- opposition_nodes %>%
  select(year, actor, system_status, parliamentary) %>%
  left_join(regime_treatment, by = c("year", "actor")) %>%
  mutate(
    access = replace_na(access, 0),
    pressure = replace_na(pressure, 0)
  ) %>%
  group_by(year) %>%
  mutate(
    n_actors = n(),
    access_mean = mean(access),
    pressure_mean = mean(pressure),
    access_sd_pop = pop_sd(access),
    pressure_sd_pop = pop_sd(pressure),
    z_access = if_else(access_sd_pop == 0 | is.na(access_sd_pop), 0, (access - access_mean) / access_sd_pop),
    z_pressure = if_else(pressure_sd_pop == 0 | is.na(pressure_sd_pop), 0, (pressure - pressure_mean) / pressure_sd_pop),
    adverse_treatment = z_pressure - z_access,
    raw_adverse_score = pressure - access,
    adverse_median = median(adverse_treatment),
    above_median_adverse = adverse_treatment > adverse_median,
    at_below_median_adverse = adverse_treatment <= adverse_median,
    raw_pressure_gt_access = pressure > access,
    raw_pressure_le_access = pressure <= access,
    any_pressure = pressure > 0,
    no_pressure = pressure == 0
  ) %>%
  ungroup()

############################
# 6. Opposition-opposition edge aggregation
############################

edge_with_status <- edges %>%
  left_join(
    node_status %>% rename(source_status = system_status, source_opp = opposition_actor),
    by = c("year", "source" = "actor")
  ) %>%
  left_join(
    node_status %>% rename(target_status = system_status, target_opp = opposition_actor),
    by = c("year", "target" = "actor")
  )

opp_opp_edges <- edge_with_status %>%
  filter(source_opp == 1, target_opp == 1) %>%
  mutate(
    actor_i = pmin(source, target),
    actor_j = pmax(source, target),
    dyad_id = paste(actor_i, actor_j, sep = " | ")
  )

negative_dyads <- opp_opp_edges %>%
  filter(sign == -1) %>%
  group_by(year, actor_i, actor_j, dyad_id) %>%
  summarise(
    negative_weight = sum(weight, na.rm = TRUE),
    negative_tie = 1L,
    .groups = "drop"
  )

positive_dyads <- opp_opp_edges %>%
  filter(sign == 1) %>%
  group_by(year, actor_i, actor_j, dyad_id) %>%
  summarise(
    positive_weight = sum(weight, na.rm = TRUE),
    positive_tie = 1L,
    .groups = "drop"
  )

# Actor-level negative conflict outcomes for H1.
actor_negative_edges <- bind_rows(
  negative_dyads %>% transmute(year, actor = actor_i, partner = actor_j, negative_weight),
  negative_dyads %>% transmute(year, actor = actor_j, partner = actor_i, negative_weight)
)

actor_negative_outcomes <- actor_negative_edges %>%
  group_by(year, actor) %>%
  summarise(
    negative_weight_sum = sum(negative_weight, na.rm = TRUE),
    unique_negative_partners = n_distinct(partner),
    .groups = "drop"
  )

actor_data <- actor_data %>%
  left_join(actor_negative_outcomes, by = c("year", "actor")) %>%
  mutate(
    negative_weight_sum = replace_na(negative_weight_sum, 0),
    unique_negative_partners = replace_na(unique_negative_partners, 0),
    conflict_centrality = negative_weight_sum / (n_actors - 1),
    unique_partner_rate = unique_negative_partners / (n_actors - 1),
    any_conflict = as.integer(unique_negative_partners > 0)
  )

############################
# 7. Dyad-cycle analysis file
############################

dyads <- opposition_nodes %>%
  group_by(year) %>%
  group_modify(make_dyads_one_year) %>%
  ungroup() %>%
  mutate(dyad_id = paste(actor_i, actor_j, sep = " | "))

actor_treatment_i <- actor_data %>%
  select(year, actor, access, pressure, z_access, z_pressure, adverse_treatment,
         raw_adverse_score, above_median_adverse, at_below_median_adverse,
         raw_pressure_gt_access, raw_pressure_le_access, any_pressure, no_pressure) %>%
  rename(
    actor_i = actor,
    access_i = access,
    pressure_i = pressure,
    z_access_i = z_access,
    z_pressure_i = z_pressure,
    adverse_treatment_i = adverse_treatment,
    raw_adverse_score_i = raw_adverse_score,
    above_median_adverse_i = above_median_adverse,
    at_below_median_adverse_i = at_below_median_adverse,
    raw_pressure_gt_access_i = raw_pressure_gt_access,
    raw_pressure_le_access_i = raw_pressure_le_access,
    any_pressure_i = any_pressure,
    no_pressure_i = no_pressure
  )

actor_treatment_j <- actor_data %>%
  select(year, actor, access, pressure, z_access, z_pressure, adverse_treatment,
         raw_adverse_score, above_median_adverse, at_below_median_adverse,
         raw_pressure_gt_access, raw_pressure_le_access, any_pressure, no_pressure) %>%
  rename(
    actor_j = actor,
    access_j = access,
    pressure_j = pressure,
    z_access_j = z_access,
    z_pressure_j = z_pressure,
    adverse_treatment_j = adverse_treatment,
    raw_adverse_score_j = raw_adverse_score,
    above_median_adverse_j = above_median_adverse,
    at_below_median_adverse_j = at_below_median_adverse,
    raw_pressure_gt_access_j = raw_pressure_gt_access,
    raw_pressure_le_access_j = raw_pressure_le_access,
    any_pressure_j = any_pressure,
    no_pressure_j = no_pressure
  )

dyad_data <- dyads %>%
  left_join(actor_treatment_i, by = c("year", "actor_i")) %>%
  left_join(actor_treatment_j, by = c("year", "actor_j")) %>%
  left_join(negative_dyads %>% select(year, dyad_id, negative_tie, negative_weight),
            by = c("year", "dyad_id")) %>%
  left_join(positive_dyads %>% select(year, dyad_id, positive_tie, positive_weight),
            by = c("year", "dyad_id")) %>%
  mutate(
    negative_tie = replace_na(negative_tie, 0L),
    positive_tie = replace_na(positive_tie, 0L),
    negative_weight = replace_na(negative_weight, 0),
    positive_weight = replace_na(positive_weight, 0),
    treatment_distance = sqrt((z_access_i - z_access_j)^2 + (z_pressure_i - z_pressure_j)^2),
    access_distance = abs(z_access_i - z_access_j),
    pressure_distance = abs(z_pressure_i - z_pressure_j),
    adverse_distance = abs(adverse_treatment_i - adverse_treatment_j),
    both_high_adverse = as.integer(adverse_treatment_i > 0 & adverse_treatment_j > 0),
    both_low_adverse = as.integer(adverse_treatment_i <= 0 & adverse_treatment_j <= 0),
    both_above_median_adverse = as.integer(above_median_adverse_i & above_median_adverse_j),
    both_at_below_median_adverse = as.integer(at_below_median_adverse_i & at_below_median_adverse_j),
    both_raw_pressure_gt_access = as.integer(raw_pressure_gt_access_i & raw_pressure_gt_access_j),
    both_raw_pressure_le_access = as.integer(raw_pressure_le_access_i & raw_pressure_le_access_j),
    both_any_pressure = as.integer(any_pressure_i & any_pressure_j),
    both_no_pressure = as.integer(no_pressure_i & no_pressure_j)
  )

############################
# 8. Lagged variables
############################

actor_lags <- actor_data %>%
  select(actor, prev_year = year,
         lag_conflict_centrality = conflict_centrality,
         lag_unique_partner_rate = unique_partner_rate,
         lag_any_conflict = any_conflict)

actor_data <- actor_data %>%
  left_join(prev_lookup, by = "year") %>%
  left_join(actor_lags, by = c("actor", "prev_year"))

actor_lag_sample <- actor_data %>%
  filter(!is.na(lag_conflict_centrality))

dyad_lags <- dyad_data %>%
  select(dyad_id, prev_year = year,
         lag_negative_tie = negative_tie,
         lag_positive_tie = positive_tie,
         lag_treatment_distance = treatment_distance,
         lag_both_high_adverse = both_high_adverse,
         lag_both_low_adverse = both_low_adverse)

dyad_data <- dyad_data %>%
  left_join(prev_lookup, by = "year") %>%
  left_join(dyad_lags, by = c("dyad_id", "prev_year")) %>%
  mutate(
    repeated_dyad = !is.na(lag_treatment_distance),
    lag_negative_tie = replace_na(lag_negative_tie, 0L),
    lag_positive_tie = replace_na(lag_positive_tie, 0L),
    lag_both_high_adverse = replace_na(lag_both_high_adverse, 0L),
    lag_both_low_adverse = replace_na(lag_both_low_adverse, 0L)
  )

dyad_repeated <- dyad_data %>%
  filter(repeated_dyad)

############################
# 9. Cross-camp E-I index
############################

ei_index <- opp_opp_edges %>%
  filter(sign == 1) %>%
  mutate(
    tie_location = if_else(source_status == target_status, "internal", "external")
  ) %>%
  group_by(year, tie_location) %>%
  summarise(total_weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = tie_location, values_from = total_weight, values_fill = 0) %>%
  mutate(
    external = if_else(is.na(external), 0, external),
    internal = if_else(is.na(internal), 0, internal),
    ei_index = (external - internal) / (external + internal)
  ) %>%
  arrange(year) %>%
  select(year, external_positive_weight = external, internal_positive_weight = internal, ei_index)

write_table(ei_index, "Table_B2_cross_camp_EI_index.csv")

############################
# 10. Save constructed analysis files
############################

write_table(actor_data, "actor_cycle_data.csv")
write_table(dyad_data, "dyad_cycle_data.csv")

############################
# 11. Appendix Section D models
############################

# Table D1: H1 actor-level models.
# Main inference rule: actor-clustered SE and actor-clustered p-value.
D1 <- bind_rows(
  run_spec("D1", "(1)", "Conflict centrality", "Main", "Adverse treatment",
           actor_data, conflict_centrality ~ adverse_treatment + factor(year),
           "adverse_treatment", "cluster", "actor"),
  run_spec("D1", "(2)", "Unique partners", "Alternative outcome", "Adverse treatment",
           actor_data, unique_partner_rate ~ adverse_treatment + factor(year),
           "adverse_treatment", "cluster", "actor"),
  run_spec("D1", "(3)", "Conflict centrality", "Raw adverse", "Raw adverse score",
           actor_data, conflict_centrality ~ raw_adverse_score + factor(year),
           "raw_adverse_score", "cluster", "actor"),
  run_spec("D1", "(4)", "Any conflict", "Binary LPM", "Adverse treatment",
           actor_data, any_conflict ~ adverse_treatment + factor(year),
           "adverse_treatment", "cluster", "actor"),
  run_spec("D1", "(5)", "Conflict centrality", "Lagged DV", "Adverse treatment",
           actor_lag_sample, conflict_centrality ~ adverse_treatment + lag_conflict_centrality + factor(year),
           "adverse_treatment", "cluster", "actor"),
  run_spec("D1", "(6)", "Unique partners", "Lagged DV", "Adverse treatment",
           actor_lag_sample, unique_partner_rate ~ adverse_treatment + lag_unique_partner_rate + factor(year),
           "adverse_treatment", "cluster", "actor")
)

# Table D2: H2 dyad-level negative-tie models.
D2 <- bind_rows(
  run_spec("D2", "(1)", "Negative tie", "Main", "Treatment distance",
           dyad_data, negative_tie ~ treatment_distance + factor(year),
           "treatment_distance", "cluster", "dyad_id"),
  run_spec("D2", "(2)", "Negative tie", "HC1 SE", "Treatment distance",
           dyad_data, negative_tie ~ treatment_distance + factor(year),
           "treatment_distance", "hc1"),
  run_spec("D2", "(3)", "Negative weight", "Weighted DV", "Treatment distance",
           dyad_data, negative_weight ~ treatment_distance + factor(year),
           "treatment_distance", "cluster", "dyad_id"),
  run_spec("D2", "(4)", "Negative tie", "Scalar distance", "Adverse distance",
           dyad_data, negative_tie ~ adverse_distance + factor(year),
           "adverse_distance", "cluster", "dyad_id"),
  run_spec("D2", "(5)", "Negative tie", "Access only", "Access distance",
           dyad_data, negative_tie ~ access_distance + factor(year),
           "access_distance", "cluster", "dyad_id"),
  run_spec("D2", "(6)", "Negative tie", "Pressure only", "Pressure distance",
           dyad_data, negative_tie ~ pressure_distance + factor(year),
           "pressure_distance", "cluster", "dyad_id")
)

# Table D3: H3 shared high adverse treatment and positive cooperation.
D3 <- bind_rows(
  run_spec("D3", "(1)", "Positive tie", "Main", "Both-high adverse",
           dyad_data, positive_tie ~ both_high_adverse + factor(year),
           "both_high_adverse", "cluster", "dyad_id"),
  run_spec("D3", "(2)", "Positive tie", "HC1 SE", "Both-high adverse",
           dyad_data, positive_tie ~ both_high_adverse + factor(year),
           "both_high_adverse", "hc1"),
  run_spec("D3", "(3)", "Positive weight", "Weighted DV", "Both-high adverse",
           dyad_data, positive_weight ~ both_high_adverse + factor(year),
           "both_high_adverse", "cluster", "dyad_id"),
  run_spec("D3", "(4)", "Positive tie", "Above median", "Both-high adverse",
           dyad_data, positive_tie ~ both_above_median_adverse + factor(year),
           "both_above_median_adverse", "cluster", "dyad_id"),
  run_spec("D3", "(5)", "Positive tie", "Raw P > A", "Both raw P > A",
           dyad_data, positive_tie ~ both_raw_pressure_gt_access + factor(year),
           "both_raw_pressure_gt_access", "cluster", "dyad_id"),
  run_spec("D3", "(6)", "Positive tie", "Any pressure", "Both any pressure",
           dyad_data, positive_tie ~ both_any_pressure + factor(year),
           "both_any_pressure", "cluster", "dyad_id"),
  permute_adverse_test("D3", "(7)", "positive_tie", "Permutation", "Both-high adverse",
                       dyad_data, actor_data,
                       observed_var = "both_high_adverse",
                       mode = "high",
                       B = 9999, seed = 42)
)

# Table D4: H4 shared low adverse treatment counter-test.
D4 <- bind_rows(
  run_spec("D4", "(1)", "Positive tie", "Main", "Both-low adverse",
           dyad_data, positive_tie ~ both_low_adverse + factor(year),
           "both_low_adverse", "cluster", "dyad_id"),
  run_spec("D4", "(2)", "Positive tie", "HC1 SE", "Both-low adverse",
           dyad_data, positive_tie ~ both_low_adverse + factor(year),
           "both_low_adverse", "hc1"),
  run_spec("D4", "(3)", "Positive weight", "Weighted DV", "Both-low adverse",
           dyad_data, positive_weight ~ both_low_adverse + factor(year),
           "both_low_adverse", "cluster", "dyad_id"),
  run_spec("D4", "(4)", "Positive tie", "At/below median", "Both at/below median",
           dyad_data, positive_tie ~ both_at_below_median_adverse + factor(year),
           "both_at_below_median_adverse", "cluster", "dyad_id"),
  run_spec("D4", "(5)", "Positive tie", "Raw P <= A", "Both raw P <= A",
           dyad_data, positive_tie ~ both_raw_pressure_le_access + factor(year),
           "both_raw_pressure_le_access", "cluster", "dyad_id"),
  run_spec("D4", "(6)", "Positive tie", "No pressure", "Both no pressure",
           dyad_data, positive_tie ~ both_no_pressure + factor(year),
           "both_no_pressure", "cluster", "dyad_id"),
  permute_adverse_test("D4", "(7)", "positive_tie", "Permutation", "Both-low adverse",
                       dyad_data, actor_data,
                       observed_var = "both_low_adverse",
                       mode = "low",
                       B = 9999, seed = 42)
)

# Table D5: dyad-level lagged controls and restricted-sample checks.
D5 <- bind_rows(
  run_spec("D5", "(1)", "Negative tie", "Lagged conflict", "Treatment distance",
           dyad_data, negative_tie ~ treatment_distance + lag_negative_tie + factor(year),
           "treatment_distance", "cluster", "dyad_id"),
  run_spec("D5", "(2)", "Negative tie", "No prior conflict", "Treatment distance",
           dyad_data %>% filter(lag_negative_tie == 0),
           negative_tie ~ treatment_distance + factor(year),
           "treatment_distance", "cluster", "dyad_id"),
  run_spec("D5", "(3)", "Positive tie", "Lagged positive tie", "Both-high adverse",
           dyad_data, positive_tie ~ both_high_adverse + lag_positive_tie + factor(year),
           "both_high_adverse", "cluster", "dyad_id"),
  run_spec("D5", "(4)", "Positive tie", "No prior positive", "Both-high adverse",
           dyad_data %>% filter(lag_positive_tie == 0),
           positive_tie ~ both_high_adverse + factor(year),
           "both_high_adverse", "cluster", "dyad_id")
)

# Table D6a: H1 reverse-direction diagnostics.
# These are reported with conventional OLS standard errors as reverse-direction diagnostics.
D6a <- bind_rows(
  run_spec("D6a", "(1)", "Current adverse treatment", "Year FE", "Lagged conflict centrality",
           actor_lag_sample, adverse_treatment ~ lag_conflict_centrality + factor(year),
           "lag_conflict_centrality", "ols"),
  run_spec("D6a", "(2)", "Current adverse treatment", "Actor FE + year FE", "Lagged conflict centrality",
           actor_lag_sample, adverse_treatment ~ lag_conflict_centrality + factor(actor) + factor(year),
           "lag_conflict_centrality", "ols")
)

# Table D6b: dyad-level reverse-direction diagnostics.
D6b <- bind_rows(
  run_spec("D6b", "(1)", "Current treatment distance", "Year FE", "Lagged negative tie",
           dyad_repeated, treatment_distance ~ lag_negative_tie + factor(year),
           "lag_negative_tie", "cluster", "dyad_id"),
  run_spec("D6b", "(2)", "Current both-high status", "Year FE", "Lagged positive tie",
           dyad_repeated, both_high_adverse ~ lag_positive_tie + factor(year),
           "lag_positive_tie", "cluster", "dyad_id")
)

############################
# 12. Descriptive comparisons used in the paper/report
############################

H1_descriptive <- actor_data %>%
  mutate(high_adverse = adverse_treatment > 0) %>%
  group_by(high_adverse) %>%
  summarise(
    n = n(),
    mean_conflict_centrality = mean(conflict_centrality),
    .groups = "drop"
  )

H2_distance_description <- dyad_data %>%
  mutate(
    top_quartile_distance = treatment_distance > quantile(treatment_distance, 0.75, na.rm = TRUE),
    at_below_median_distance = treatment_distance <= median(treatment_distance, na.rm = TRUE)
  ) %>%
  summarise(
    top_quartile_n = sum(top_quartile_distance),
    top_quartile_negative_rate = mean(negative_tie[top_quartile_distance]),
    at_below_median_n = sum(at_below_median_distance),
    at_below_median_negative_rate = mean(negative_tie[at_below_median_distance]),
    relative_ratio = top_quartile_negative_rate / at_below_median_negative_rate
  )

H3_descriptive <- dyad_data %>%
  group_by(both_high_adverse) %>%
  summarise(
    dyad_cycle_n = n(),
    positive_tie_dyads = sum(positive_tie),
    positive_tie_rate = mean(positive_tie),
    .groups = "drop"
  )

H4_descriptive <- dyad_data %>%
  group_by(both_low_adverse) %>%
  summarise(
    dyad_cycle_n = n(),
    positive_tie_dyads = sum(positive_tie),
    positive_tie_rate = mean(positive_tie),
    .groups = "drop"
  )

write_table(H1_descriptive, "H1_descriptive_high_low_adverse.csv")
write_table(H2_distance_description, "H2_descriptive_distance_comparison.csv")
write_table(H3_descriptive, "H3_descriptive_both_high.csv")
write_table(H4_descriptive, "H4_descriptive_both_low.csv")

############################
# 13. Export Appendix D results
############################

appendix_D_results <- bind_rows(D1, D2, D3, D4, D5, D6a, D6b)

write_table(D1, "Table_D1_H1_actor_level_models.csv")
write_table(D2, "Table_D2_H2_negative_tie_models.csv")
write_table(D3, "Table_D3_H3_high_adverse_positive_cooperation.csv")
write_table(D4, "Table_D4_H4_low_adverse_counter_test.csv")
write_table(D5, "Table_D5_dyad_lagged_restricted_checks.csv")
write_table(D6a, "Table_D6a_H1_reverse_direction.csv")
write_table(D6b, "Table_D6b_dyad_reverse_direction.csv")
write_table(appendix_D_results, "Appendix_D_all_results_long.csv")

# A compact display table with formatted coefficient / SE / p-value.
appendix_D_display <- appendix_D_results %>%
  select(table, model, outcome, specification, predictor_label,
         coef_display, se_display, p_display, n, inference_rule, year_fe)
write_table(appendix_D_display, "Appendix_D_display_table_long.csv")

############################
# 14. Console summary
############################

cat("\nReplication complete.\n")
cat("Output folder:", normalizePath(output_dir), "\n\n")
cat("Key constructed files:\n")
cat("  - actor_cycle_data.csv\n")
cat("  - dyad_cycle_data.csv\n")
cat("  - Appendix_D_all_results_long.csv\n")
cat("  - Appendix_D_display_table_long.csv\n\n")

cat("Dataset scope check:\n")
print(scope_by_cycle)

cat("\nAppendix D display table:\n")
print(appendix_D_display, n = nrow(appendix_D_display))

############################################################
# End of script
############################################################
