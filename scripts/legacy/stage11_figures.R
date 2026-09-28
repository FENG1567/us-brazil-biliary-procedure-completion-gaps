#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(svglite); library(ragg); library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sourcedir <- file.path(outdir, "source_data"); dir.create(sourcedir, recursive = TRUE, showWarnings = FALSE)
res <- file.path(root, "results", "stage11"); qc <- file.path(root, "qc", "stage11")

blue <- "#2C6EAA"; teal <- "#1B8A7A"; orange <- "#D9772A"; red <- "#B84A4A"
grey <- "#6B7280"; light <- "#E5E7EB"; dark <- "#202124"; purple <- "#7561A8"
cohort_labs <- c(acute_cholecystitis = "Acute cholecystitis", choledocholithiasis = "Choledocholithiasis")
trt_cols <- c("Early completion" = blue, "Not completed by day 3" = grey)

theme_pub <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = .32, colour = dark),
          axis.ticks = element_line(linewidth = .32, colour = dark),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - .3, colour = dark),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = base_size),
          legend.title = element_blank(), legend.text = element_text(size = base_size - .4),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - .2, colour = grey),
          plot.margin = margin(4, 5, 4, 5), panel.grid = element_blank())
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 178) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite(file.path(outdir, paste0(stem, ".svg")), width = w, height = h, system_fonts = list(Arial = "Arial")); print(plot); dev.off()
  cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = w, height = h, family = "Arial"); print(plot); dev.off()
  agg_tiff(file.path(outdir, paste0(stem, ".tiff")), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(plot); dev.off()
  agg_png(file.path(outdir, paste0(stem, ".png")), width = w, height = h, units = "in", res = 600, background = "white"); print(plot); dev.off()
}

tags <- plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", family = "Arial", size = 8)))

# Figure 1: evidence architecture, attrition, fixed time zero, and adjustment logic.
cap <- CJ(capability = c("National weighting", "Within-year linkage", "Inpatient cost", "Structural resources", "Geographic care flow", "Clinical physiology"),
          database = c("NRD", "SIH-SUS/CNES", "MIMIC-IV"))
cap[, supported := (database == "NRD" & capability %chin% c("National weighting", "Within-year linkage", "Inpatient cost")) |
      (database == "SIH-SUS/CNES" & capability %chin% c("Structural resources", "Geographic care flow")) |
      (database == "MIMIC-IV" & capability %chin% c("Within-year linkage", "Clinical physiology"))]
cap[, database := factor(database, levels = c("NRD", "SIH-SUS/CNES", "MIMIC-IV"))]
cap[, capability := factor(capability, levels = rev(c("National weighting", "Within-year linkage", "Inpatient cost", "Structural resources", "Geographic care flow", "Clinical physiology")))]
fwrite(cap, file.path(sourcedir, "Figure1A_database_capabilities.csv"))
p1a <- ggplot(cap, aes(database, capability)) + geom_tile(fill = "white", color = light, linewidth = .35) +
  geom_point(aes(fill = supported), shape = 21, size = 3.0, stroke = .5) +
  scale_fill_manual(values = c(`TRUE` = blue, `FALSE` = "white")) +
  labs(title = "Complementary evidence roles", subtitle = "Patient-level records were not pooled", x = NULL, y = NULL) +
  theme_pub() + theme(legend.position = "none", axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(face = "bold"))

flow <- fread(file.path(qc, "nrd_day3_flow.csv"))
flow_long <- melt(flow, id.vars = "cohort", measure.vars = c("source_records", "principal_diagnosis", "principal_urgent", "no_contraindication_proxy", "primary_day3_eligible"),
                  variable.name = "stage", value.name = "n")
flow_long[, stage := factor(stage, levels = c("source_records", "principal_diagnosis", "principal_urgent", "no_contraindication_proxy", "primary_day3_eligible"),
                          labels = c("Index cohort", "Principal diagnosis", "Urgent", "No coded proxy", "Day-3 eligible"))]
flow_long[, cohort_label := cohort_labs[cohort]]
fwrite(flow_long, file.path(sourcedir, "Figure1B_nrd_flow.csv"))
p1b <- ggplot(flow_long, aes(stage, n, color = cohort_label, group = cohort_label)) +
  geom_line(linewidth = .65) + geom_point(size = 2.0) +
  geom_text(aes(label = comma(n)), vjust = -0.65, size = 1.85, family = "Arial", show.legend = FALSE) +
  scale_color_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), expand = expansion(mult = c(.02, .15))) +
  labs(title = "NRD analytic attrition", x = NULL, y = "Admissions") + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

timeline <- data.table(day = c(0, 3, 4, 93), pos = c(0, 1.6, 3.3, 10),
                       y = c(.82, .20, .82, .82),
                       label = c("Admission", "Day 3 landmark\nclassify early treatment", "Day 4\nfollow-up starts", "Day 93\nfollow-up ends"))
fwrite(timeline, file.path(sourcedir, "Figure1C_landmark_timeline.csv"))
p1c <- ggplot() +
  annotate("segment", x = 0, xend = 10, y = .55, yend = .55, linewidth = .75, colour = dark, arrow = arrow(length = unit(.10, "inches"))) +
  annotate("rect", xmin = 0, xmax = 1.6, ymin = .43, ymax = .67, fill = alpha(orange, .28), colour = orange, linewidth = .35) +
  annotate("rect", xmin = 3.3, xmax = 10, ymin = .43, ymax = .67, fill = alpha(blue, .16), colour = blue, linewidth = .35) +
  geom_point(data = timeline, aes(pos, .55), size = 2.0, color = c(dark, orange, blue, blue)) +
  geom_text(data = timeline, aes(pos, y, label = label), size = 2.15, family = "Arial", lineheight = .92) +
  annotate("text", x = 6.8, y = .18, label = "Biliary readmission/death outcomes;\nlater treatment does not reassign exposure", size = 1.9, family = "Arial") +
  coord_cartesian(xlim = c(-.35, 10.35), ylim = c(.05, 1.00), clip = "off") +
  labs(title = "Fixed day-3 landmark and common time zero") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 7.5), plot.margin = margin(7, 10, 5, 10))

nodes <- data.table(x = c(1, 1, 3, 3, 5), y = c(3, 1, 3, 1, 2),
                    label = c("Baseline disease and\ncandidacy proxies", "Hospital and access\ncharacteristics", "Completion by\nday 3", "Post-treatment course\nand length of stay", "Day 4-93\noutcomes"),
                    type = c("adjust", "adjust", "exposure", "post", "outcome"))
edges <- data.table(x = c(1.55, 1.55, 3.55, 3.55, 1.55), y = c(3, 1, 3, 1, 3),
                    xend = c(2.45, 2.45, 4.45, 4.45, 4.45), yend = c(3, 3, 2, 2, 2))
fwrite(nodes, file.path(sourcedir, "Figure1D_adjustment_nodes.csv")); fwrite(edges, file.path(sourcedir, "Figure1D_adjustment_edges.csv"))
p1d <- ggplot() + geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(.08, "inches")), linewidth = .45, color = grey) +
  geom_label(data = nodes, aes(x, y, label = label, fill = type), size = 2.1, family = "Arial", label.size = .3, label.padding = unit(.10, "lines")) +
  scale_fill_manual(values = c(adjust = "#EAF2F8", exposure = "#FCE8D5", post = "#F3F4F6", outcome = "#E6F3F0")) +
  annotate("text", x = 3, y = .30, label = "Only admission-recorded variables entered adjustment;\npost-treatment variables were not adjusted", size = 1.95, family = "Arial", color = red) +
  coord_cartesian(xlim = c(.25, 5.75), ylim = c(0, 3.65), clip = "off") + labs(title = "Prespecified adjustment boundary") + theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 7.5), plot.margin = margin(6, 8, 5, 8))

fig1 <- (p1a | p1b) / (p1c | p1d) + plot_layout(heights = c(1.08, .92)) + tags
save_pub(fig1, "Figure1_Study_Architecture")

# Figure 2: primary NRD associations, bias signal, balance, sensitivity.
primary <- fread(file.path(res, "nrd_day3_primary_estimates.csv"))
boot <- fread(file.path(res, "nrd_day3_design_bootstrap_summary.csv"))
primary[, cohort_label := cohort_labs[cohort]]
risk <- melt(primary[outcome %in% c("biliary", "composite")],
             id.vars = c("cohort", "cohort_label", "outcome"), measure.vars = c("risk_early", "risk_not_early"),
             variable.name = "treatment_group", value.name = "risk")
risk[, treatment_group := factor(treatment_group, levels = c("risk_not_early", "risk_early"), labels = c("Not completed by day 3", "Early completion"))]
risk[, outcome_label := fifelse(outcome == "biliary", "Biliary readmission", "Biliary readmission or inpatient death")]
fwrite(risk, file.path(sourcedir, "Figure2A_weighted_risks.csv"))
p2a <- ggplot(risk, aes(treatment_group, risk, color = treatment_group, group = outcome_label)) +
  geom_line(aes(group = interaction(cohort_label, outcome_label)), color = light, linewidth = .6) + geom_point(size = 2.2) +
  facet_grid(outcome_label ~ cohort_label) + scale_color_manual(values = trt_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(.02, .12))) +
  labs(title = "Overlap-weighted absolute risks", x = NULL, y = "Day 4-93 risk") + theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

rd <- primary[, .(cohort, cohort_label, outcome, estimate = risk_difference, lcl = rd_lcl, ucl = rd_ucl)]
rd[, outcome_label := factor(outcome, levels = c("biliary", "composite", "death", "nonbiliary"),
                             labels = c("Biliary readmission", "Biliary or inpatient death", "Observed inpatient death", "Nonbiliary readmission"))]
fwrite(rd, file.path(sourcedir, "Figure2B_risk_differences.csv"))
p2b <- ggplot(rd, aes(estimate, outcome_label, color = outcome_label)) + geom_vline(xintercept = 0, color = grey, linewidth = .35) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = .15, linewidth = .55) + geom_point(size = 2.1) +
  facet_wrap(~cohort_label, scales = "free_x") +
  scale_color_manual(values = c("Biliary readmission" = red, "Biliary or inpatient death" = orange, "Observed inpatient death" = purple, "Nonbiliary readmission" = grey)) +
  scale_x_continuous(labels = percent_format(accuracy = .1)) +
  labs(title = "Adjusted risk differences", subtitle = "Survey-linearized 95% confidence intervals", x = "Early minus not-early completion", y = NULL) + theme_pub() + theme(legend.position = "none")

bal <- fread(file.path(qc, "nrd_day3_balance.csv")); bal[, cohort_label := cohort_labs[cohort]]
top_vars <- bal[, .(max_pre = max(abs(smd_before))), by = variable][order(-max_pre)][1:min(16, .N), variable]
bal <- bal[variable %chin% top_vars]
bal_long <- melt(bal, id.vars = c("cohort", "cohort_label", "variable"), measure.vars = c("smd_before", "smd_after"), variable.name = "weighting", value.name = "smd")
bal_long[, weighting := factor(weighting, levels = c("smd_before", "smd_after"), labels = c("Before weighting", "After overlap weighting"))]
fwrite(bal_long, file.path(sourcedir, "Figure2C_balance.csv"))
p2c <- ggplot(bal_long, aes(abs(smd), reorder(variable, abs(smd)), color = weighting)) + geom_vline(xintercept = .10, linetype = 2, color = red, linewidth = .35) +
  geom_point(size = 1.4) + facet_wrap(~cohort_label) + scale_color_manual(values = c("Before weighting" = grey, "After overlap weighting" = blue)) +
  scale_x_continuous(limits = c(0, max(.52, max(abs(bal_long$smd)) * 1.03))) +
  labs(title = "Measured-covariate balance", x = "Absolute standardized mean difference", y = NULL) + theme_pub(6.1) + theme(legend.position = "bottom")

sens <- fread(file.path(res, "nrd_day3_sensitivity_estimates.csv")); sens[, cohort_label := cohort_labs[cohort]]
sens[, variant_label := factor(variant, levels = c("broad_eligibility_day3", "age18_79_day3", "exclude_transfer_day3", "principal_all_admissions_day3", "primary_day2", "primary_day5", "primary_day7", "strict_procedure_day3", "broad_procedure_day3"),
                               labels = c("Include coded proxies", "Age 18-79", "Exclude transfers", "All admission routes", "Day-2 landmark", "Day-5 landmark", "Day-7 landmark", "Total chole only", "Broad therapeutic ERCP"))]
fwrite(sens, file.path(sourcedir, "Figure2D_sensitivity.csv"))
p2d <- ggplot(sens, aes(risk_difference, variant_label, color = cohort_label)) + geom_vline(xintercept = 0, color = grey, linewidth = .35) +
  geom_errorbarh(aes(xmin = rd_lcl, xmax = rd_ucl), height = .12, linewidth = .45,
                 position = position_dodge(width = .45)) +
  geom_point(size = 1.7, position = position_dodge(width = .45)) +
  scale_color_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  scale_x_continuous(labels = percent_format(accuracy = .1), n.breaks = 4) +
  labs(title = "Prespecified design sensitivities", x = "Biliary readmission risk difference", y = NULL) +
  theme_pub(6.1) + theme(legend.position = "bottom")

fig2 <- (p2a | p2b) / (p2c | p2d) + plot_layout(heights = c(.92, 1.25)) + tags
save_pub(fig2, "Figure2_Landmark_Associations_and_Bias")

# Figure 3: recurrence phenotype and aligned scenarios.
sub <- fread(file.path(res, "nrd_day3_subtype_estimates.csv")); sub[, cohort_label := cohort_labs[cohort]]
fwrite(sub, file.path(sourcedir, "Figure3A_recurrence_subtypes.csv"))
p3a <- ggplot(sub, aes(risk_difference, reorder(outcome, risk_difference), color = cohort_label)) + geom_vline(xintercept = 0, color = grey, linewidth = .35) +
  geom_errorbarh(aes(xmin = rd_lcl, xmax = rd_ucl), height = .14, linewidth = .5) + geom_point(size = 1.9) +
  facet_wrap(~cohort_label, scales = "free_y") + scale_color_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  scale_x_continuous(labels = percent_format(accuracy = .1)) + labs(title = "Recurrence phenotypes", x = "Risk difference", y = NULL) + theme_pub() + theme(legend.position = "none")

bur <- fread(file.path(res, "nrd_day3_aligned_burden_point.csv")); bur[, cohort_label := cohort_labs[cohort]]
bur[, closure := factor(scenario, levels = c("25_percent_of_supported_not_early_group", "50_percent_of_supported_not_early_group", "75_percent_of_supported_not_early_group"), labels = c("25%", "50%", "75%"))]
bur[, pct := as.integer(sub("_percent.*", "", scenario))]
bur[, `:=`(burden_metric = paste0("burden", pct), cost_metric = paste0("netcost", pct))]
bur_ci <- merge(
  boot[grepl("^burden", metric), .(cohort, burden_metric = metric, burden_lcl = lcl, burden_ucl = ucl)],
  boot[grepl("^netcost", metric), .(cohort, cost_metric = metric, cost_lcl = lcl, cost_ucl = ucl)],
  by = "cohort", allow.cartesian = TRUE
)
bur <- merge(bur, bur_ci, by = c("cohort", "burden_metric", "cost_metric"), all.x = TRUE)
fwrite(bur, file.path(sourcedir, "Figure3B_C_aligned_burden.csv"))
p3b <- ggplot(bur, aes(closure, associated_biliary_readmissions, fill = cohort_label)) + geom_col(position = position_dodge(width = .72), width = .62) +
  geom_errorbar(aes(ymin = burden_lcl, ymax = burden_ucl), position = position_dodge(width = .72), width = .16, linewidth = .4) +
  geom_text(aes(label = comma(round(associated_biliary_readmissions))), position = position_dodge(width = .72), vjust = -.25, size = 1.9, family = "Arial") +
  scale_fill_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .14))) +
  labs(title = "Associated biliary admissions", subtitle = "Same eligible untreated common-support population", x = "Scenario uptake", y = "Admissions") + theme_pub() + theme(legend.position = "bottom")

p3c <- ggplot(bur, aes(closure, net_inpatient_facility_cost_change_2022usd / 1e6, fill = cohort_label)) + geom_hline(yintercept = 0, color = grey, linewidth = .35) +
  geom_col(position = position_dodge(width = .72), width = .62) +
  geom_errorbar(aes(ymin = cost_lcl / 1e6, ymax = cost_ucl / 1e6), position = position_dodge(width = .72), width = .16, linewidth = .4) +
  geom_text(aes(label = dollar(net_inpatient_facility_cost_change_2022usd / 1e6, accuracy = .1, suffix = "M")), position = position_dodge(width = .72),
            vjust = ifelse(bur$net_inpatient_facility_cost_change_2022usd >= 0, -.25, 1.2), size = 1.8, family = "Arial") +
  scale_fill_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  labs(title = "Net 93-day inpatient facility cost", subtitle = "Index plus biliary readmission costs; 2022 US dollars", x = "Scenario uptake", y = "Change, millions") + theme_pub() + theme(legend.position = "bottom")

qba <- fread(file.path(res, "nrd_day3_quantitative_bias_grid.csv")); qba <- qba[prevalence_early == .10 & prevalence_not_early %in% c(.20, .40, .60)]
qba[, cohort_label := cohort_labs[cohort]]; qba[, prevalence_gap := factor(prevalence_not_early - prevalence_early)]
fwrite(qba, file.path(sourcedir, "Figure3D_unmeasured_confounding.csv"))
p3d <- ggplot(qba, aes(factor(confounder_outcome_rr), prevalence_gap, fill = corrected_rr)) + geom_tile(color = "white", linewidth = .35) +
  geom_text(aes(label = sprintf("%.2f", corrected_rr)), size = 2.0, family = "Arial") + facet_wrap(~cohort_label) +
  scale_fill_gradient2(low = teal, mid = "white", high = red, midpoint = 1,
                       limits = c(min(1, qba$corrected_rr), max(1, qba$corrected_rr))) +
  labs(title = "Unmeasured-candidacy sensitivity", subtitle = "Corrected early/not-early risk ratio", x = "Outcome RR for unmeasured factor", y = "Prevalence excess in not-early group") + theme_pub(6.2) + theme(legend.position = "bottom")

fig3 <- p3a / (p3b | p3c) / p3d + plot_layout(heights = c(1.05, .92, .93)) + tags
save_pub(fig3, "Figure3_Aligned_Burden_and_Recurrence", 183, 210)

# Figure 4: Brazil care flow, two-way clustered access associations, coding bounds, MIMIC ancillary analysis.
bflow <- fread(file.path(res, "brazil_care_flow_descriptives.csv")); bflow[, cohort_label := cohort_labs[cohort]]
bflow[, share := n / sum(n), by = cohort]; bflow[, flow_label := factor(care_flow_category, levels = c("same_municipality", "within_state_external", "interstate"), labels = c("Same municipality", "Within-state external", "Interstate"))]
fwrite(bflow, file.path(sourcedir, "Figure4A_brazil_care_flow.csv"))
p4a <- ggplot(bflow, aes(flow_label, share, fill = flow_label)) + geom_col(width = .64) +
  geom_text(aes(label = percent(share, accuracy = .1)), vjust = -.25, size = 1.9, family = "Arial") + facet_wrap(~cohort_label) +
  scale_fill_manual(values = c("Same municipality" = light, "Within-state external" = teal, "Interstate" = orange)) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, .13))) + labs(title = "Observed residence-to-hospital care flow", x = NULL, y = "Admissions") + theme_pub() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")

bcon <- fread(file.path(res, "brazil_access_marginal_contrasts.csv")); bcon[, cohort_label := cohort_labs[cohort]]
bcon[, contrast_label := fifelse(variable == "structural_readiness", "Structural readiness contrast",
                          fifelse(variable == "local_prior_capability", "Residence-area prior capability",
                          fifelse(variable == "log_destination_prior_operations", "Destination prior operations contrast",
                          fifelse(high == "within_state_external", "Within-state external vs local", "Interstate vs local"))))]
fwrite(bcon, file.path(sourcedir, "Figure4B_brazil_contrasts.csv"))
p4b <- ggplot(bcon, aes(risk_difference, reorder(contrast_label, risk_difference), color = cohort_label)) + geom_vline(xintercept = 0, color = grey, linewidth = .35) +
  geom_errorbarh(aes(xmin = rd_lcl, xmax = rd_ucl), height = .14, linewidth = .5,
                 position = position_dodge(width = .48)) +
  geom_point(size = 1.9, position = position_dodge(width = .48)) +
  scale_color_manual(values = c("Acute cholecystitis" = blue, "Choledocholithiasis" = teal)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), n.breaks = 5) +
  labs(title = "Brazilian access associations", subtitle = "Residence and destination two-way clustered 95% CIs", x = "Completion risk difference", y = NULL) +
  theme_pub(6.2) + theme(legend.position = "bottom")

cross <- fread(file.path(res, "cross_system_2021_2022_coding_bounds.csv")); cross[, cohort_label := cohort_labs[cohort]]
cross_rate <- dcast(cross[estimand %in% c("us_lower", "us_upper", "br_lower", "br_upper")], cohort + cohort_label ~ estimand, value.var = "estimate")
cross_long <- rbindlist(list(
  cross_rate[, .(cohort, cohort_label, system = "United States", lower = us_lower, upper = us_upper)],
  cross_rate[, .(cohort, cohort_label, system = "Brazil", lower = br_lower, upper = br_upper)]
))
fwrite(cross_long, file.path(sourcedir, "Figure4C_cross_system_bounds.csv"))
p4c <- ggplot(cross_long, aes(y = system, color = system)) + geom_segment(aes(x = lower, xend = upper, yend = system), linewidth = 2.4, alpha = .55) +
  geom_point(aes(x = lower), shape = 21, fill = "white", size = 2.0) + geom_point(aes(x = upper), shape = 16, size = 2.0) +
  facet_wrap(~cohort_label) + scale_color_manual(values = c("United States" = blue, "Brazil" = teal)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) + labs(title = "2021-2022 coding-bound completion", subtitle = "Open/filled points: lower/upper code definitions", x = "Overlap-population completion", y = NULL) + theme_pub() + theme(legend.position = "none")

mimic <- fread(file.path(res, "mimic_day3_landmark_estimates.csv"))
mimic <- mimic[outcome == "biliary" & stratum == "all"]
mimic[, cohort_label := cohort_labs[cohort]]
mimic[, model_label := factor(model, levels = c("admission_recorded_covariates", "earliest_recorded_labs_sensitivity"), labels = c("Admission-recorded", "+ earliest labs"))]
fwrite(mimic, file.path(sourcedir, "Figure4D_mimic_ancillary.csv"))
p4d <- ggplot(mimic, aes(risk_difference, model_label, color = model_label)) + geom_vline(xintercept = 0, color = grey, linewidth = .35) +
  geom_errorbarh(aes(xmin = rd_lcl, xmax = rd_ucl), height = .14, linewidth = .5) + geom_point(size = 1.9) + facet_wrap(~cohort_label, scales = "free_x") +
  scale_color_manual(values = c("Admission-recorded" = grey, "+ earliest labs" = blue)) + scale_x_continuous(labels = percent_format(accuracy = .1)) +
  labs(title = "MIMIC-IV ancillary landmark analysis", subtitle = "Single system; laboratory timing may follow treatment", x = "Biliary readmission risk difference", y = NULL) + theme_pub(6.2) + theme(legend.position = "bottom")

fig4 <- (p4a | p4b) / (p4c | p4d) + plot_layout(heights = c(1.12, .88)) + tags
save_pub(fig4, "Figure4_Care_Flow_and_Cross_System_Bounds")

cat("FIGURES_PASS four_pdf_four_svg_four_tiff_four_png\n")
