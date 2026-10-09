library(SSN2)
library(ggplot2)


bugs <- ssn_import("bugs.ssn", predpts = "pred")

ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = bugs$preds$pred, color = "black", size = 0.9) +
  geom_sf(data = bugs$obs, aes(color = Rich), size = 2) +
  scale_color_viridis_c(option = "H", limits = c(0, 59)) +
  theme_bw(base_size = 16)

# ssn_create_distmat(bugs, predpts = "pred", among_predpts = TRUE)

standard_start <- Sys.time()
ssn_mod <- ssn_lm(
  formula = Rich ~ ELEV_DEM,
  ssn.object = bugs,
  tailup_type = "exponential",
  taildown_type = "exponential",
  additive = "afvArea"
)
standard_time <- Sys.time() - standard_start
standard_time

aug_ssn_mod <- augment(ssn_mod, newdata = "pred", interval = "prediction")

# ssn_create_bigdist(bugs, predpts = "pred", among_predpts = TRUE, no_cores = 2)

set.seed(1)
bd_start <- Sys.time()
ssn_mod_bd <- ssn_lm(
  formula = Rich ~ ELEV_DEM,
  ssn.object = bugs,
  tailup_type = "exponential",
  taildown_type = "exponential",
  additive = "afvArea",
  local = TRUE
)
bd_time <- Sys.time() - bd_start
bd_time

bugs$obs$index <- as.factor(ssn_mod_bd$local_index)
ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = bugs$obs, aes(color = index), size = 2) +
  scale_color_viridis_d() +
  theme_bw(base_size = 16)

cov_pred_by_obs <- covmatrix(ssn_mod_bd, newdata = "pred")
dim(cov_pred_by_obs)

avg_cov_with_obs <- colMeans(cov_pred_by_obs)

cov_order <- order(avg_cov_with_obs, decreasing = TRUE)
largest_300 <- cov_order[1:250]
row_number <- seq(from = 1, to = NROW(bugs$obs))
bugs$obs$Include <- ifelse(row_number %in% largest_300, "Yes", "No")

aug_ssn_mod_bd <- augment(
  x = ssn_mod_bd,
  newdata = "pred",
  interval = "prediction",
  local = list(size = 250)
)

summary(ssn_mod)
summary(ssn_mod_bd)

loocv(ssn_mod)
loocv(ssn_mod_bd)

head(aug_ssn_mod[, c("ELEV_DEM", ".fitted", ".lower", ".upper", "pid")])
head(aug_ssn_mod_bd[, c("ELEV_DEM", ".fitted", ".lower", ".upper", "pid")])

aug_ssn_mod$type <- "Standard"
aug_ssn_mod_bd$type <- "LNBH"
augs <- rbind(aug_ssn_mod, aug_ssn_mod_bd)
augs$type <- factor(augs$type, levels = c("Standard", "LNBH"))
ggplot(augs) +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = augs, aes(color = .fitted), size = 2) +
  facet_wrap(~ type) +
  scale_color_viridis_c(option = "H", limits = c(0, 59)) +
  scale_x_continuous(breaks = -c(117, 114)) +
  theme_bw(base_size = 16)

predict(
  object = ssn_mod,
  newdata = "pred",
  block = TRUE,
  interval = "prediction"
)
predict(
  object = ssn_mod_bd,
  newdata = "pred",
  block = TRUE,
  interval = "prediction",
  local = list(size = 250)
)
