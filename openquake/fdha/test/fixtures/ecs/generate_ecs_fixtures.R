#!/usr/bin/env Rscript
# =============================================================================
# DEV-ONLY validation oracle (NOT shipped, NOT imported at runtime).  v2.
#
# v1 dumped lpmatrix WITHOUT the u-values it was evaluated at -> unusable for
# cell-by-cell basis validation, and omitted the selection-gate and GC2-gate
# fixtures. v2 fixes the anchor and emits everything so a SINGLE R run closes
# the selection gate, the spline gate, AND the GC2 (+50 km extension) gate.
#
# Assembly block is verbatim from CalcEventCoordinateSystem.R (L37-202),
# Lavrentiadis et al. (2024) ECS reference impl. Only the emit section is new.
#
# Emits (all under <out_dir>):
#   counts.csv          n_disp,n_rup,n_data4ecs,n_replicated,min_wt  [row-count gate]
#   data4ecs_pre.csv    post-filter, pre-replication pts + weights   [selection/weight gate]
#   ecs_trace.csv       final ECS trace (ecs_main out[[1]])          [GC2 ref line + end-to-end]
#   fault_disp.csv      GAM model frame: Longitude,Latitude,wt,x,y,u,t
#                         -> u anchors lpmatrix; (lon,lat)->(u,t) is GC2 truth
#   lpmatrix.csv        mgcv tp basis at fault_disp$u                [spline basis]
#   coef.csv            fitted coefficients, named                   [spline solution]
#   penalty_S_pred1.csv / penalty_S_pred2.csv  penalty matrices + sp [spline penalty]
#
# GC2 gate usage (Python side): build MultiLine from ecs_trace[,(Longitude,
# Latitude)] as the reference line; transform fault_disp[,(Longitude,Latitude)];
# compare the resulting (u,t) to fault_disp$u / fault_disp$t (R, gc2ext_ut w/
# 50 km extension). Assert SIGN and ORIGIN, not just magnitude.
#
# Usage:
#   Rscript generate_ecs_fixtures.R \
#       <ecs_functions.R> <FLATFILE_MEASUREMENTS.csv> <FLATFILE_RUPTURES.csv> <out_dir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript generate_ecs_fixtures.R <ecs_functions.R> <measurements.csv> <ruptures.csv> <out_dir>")
}
fname_funcs <- args[1]; fname_disp <- args[2]; fname_rup <- args[3]; out_dir <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(plyr)
library(assertthat)
source(fname_funcs)   # ecs_main, weight_ecs_data, rup_length, ...
# override longlat2UTM / UTM2longlat with a self-contained TM (avoids sf);
# validated against pyproj to < 1e-3 m. Sourced AFTER ecs_functions.R.
source(file.path(dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "utm_shim.R"))

# --- parameters, verbatim from CalcEventCoordinateSystem.R (L37-54) -----------
rank2consider         <- c('Total', 'Principal', 'Cumulative')
lambda_p              <- 0.05
wt_disp_min           <- 0.05
flag_wt_rup_opt       <- 2
abs_wt_rup            <- 0.05
flag_wt_rup_len       <- TRUE
flag_norm_disp_rup_wt <- FALSE
r_thres_max           <- 1e4
field_disp_wt         <- 'recommended_net_preferred_for_analysis_meters'

# --- read flatfiles, normalise lon/lat column names (L62-71) ------------------
flatfile_disp <- read.csv(fname_disp, header = TRUE, sep = ",")
flatfile_rup  <- read.csv(fname_rup,  header = TRUE, sep = ",")
colnames(flatfile_rup)[colnames(flatfile_rup)  == "longitude_degrees"] <- "Longitude"
colnames(flatfile_rup)[colnames(flatfile_rup)  == "latitude_degrees"]  <- "Latitude"
colnames(flatfile_disp)[colnames(flatfile_disp) == "longitude_degrees"] <- "Longitude"
colnames(flatfile_disp)[colnames(flatfile_disp) == "latitude_degrees"]  <- "Latitude"

# single-event fixture: FAIL LOUDLY if the flatfile holds more than one event
assert_that(length(unique(flatfile_disp$EQ_ID)) == 1,
            msg = "flatfile has >1 EQ_ID; this oracle is single-event. Subset first.")
eqid <- unique(flatfile_disp$EQ_ID)[1]
message("Event EQ_ID: ", eqid)
fault_disp_all <- subset(flatfile_disp, EQ_ID == eqid)
fault_rup_all  <- subset(flatfile_rup,  EQ_ID == eqid)

# --- outlier removal (L120-124) ----------------------------------------------
i_outliers_d <- fault_disp_all[, field_disp_wt] < -900
fault_disp <- fault_disp_all[!i_outliers_d, ]
fault_rup  <- fault_rup_all

# --- displacement weights (L128-139) -----------------------------------------
wt_array_disp <- abs(fault_disp[, field_disp_wt])
i_disp2keep   <- (!is.na(wt_array_disp)) & (fault_disp$rank %in% rank2consider)
wt_array_disp <- wt_array_disp[i_disp2keep]
wt_array_disp[wt_array_disp < wt_disp_min] <- wt_disp_min
min_wt_array_disp <- min(wt_array_disp, na.rm = TRUE)

# --- rupture weights (L140-163, opt 2 + length scaling) ----------------------
wt_array_rup    <- array(abs_wt_rup, dim = nrow(fault_rup))
wt_array_ruplen <- array(dim = nrow(fault_rup))
for (r_id in unique(fault_rup$RUP_ID)) {
  i_rup   <- fault_rup$RUP_ID == r_id
  seg_rup <- fault_rup[i_rup, c('Longitude', 'Latitude', 'RUP_ID', 'NODE_ID')]
  seg_len <- rup_length(seg_rup)
  wt_array_ruplen[i_rup] <- seg_len / sum(i_rup)
}
if (flag_wt_rup_len) wt_array_rup <- wt_array_rup * wt_array_ruplen
i_rup2keep   <- (!is.na(wt_array_rup)) & (fault_rup$rank %in% rank2consider)
wt_array_rup <- wt_array_rup[i_rup2keep]

# --- normalise weights (L165-177) --------------------------------------------
min_wt        <- min(min(wt_array_rup), min(wt_array_disp))
wt_array_disp <- round(wt_array_disp / min_wt)
wt_array_rup  <- round(wt_array_rup  / min_wt)

data4ecs1 <- fault_disp[i_disp2keep, c('Longitude', 'Latitude')]; data4ecs1$wt <- wt_array_disp
data4ecs2 <- fault_rup[i_rup2keep, c('Longitude', 'Latitude', 'RUP_ID', 'NODE_ID')]; data4ecs2$wt <- wt_array_rup
data4ecs    <- rbind.fill(data4ecs1, data4ecs2)
data4ecs_wt <- weight_ecs_data(data4ecs, ratio_thres = r_thres_max)

# =============================================================================
# EMIT (v2 additions)
# =============================================================================

# --- (1) selection / row-count gate ------------------------------------------
write.csv(data4ecs, file.path(out_dir, "data4ecs_pre.csv"), row.names = FALSE)
write.csv(data.frame(n_disp_kept  = nrow(data4ecs1),
                     n_rup_kept   = nrow(data4ecs2),
                     n_data4ecs   = nrow(data4ecs),
                     n_replicated = nrow(data4ecs_wt),
                     min_wt       = min_wt),
          file.path(out_dir, "counts.csv"), row.names = FALSE)
message("counts: disp=", nrow(data4ecs1), " rup=", nrow(data4ecs2),
        " data4ecs=", nrow(data4ecs), " replicated=", nrow(data4ecs_wt))

# --- compute ECS (L202) ------------------------------------------------------
out        <- ecs_main(data4ecs_wt, lamb_spline = lambda_p, calc_xy = TRUE)
ecs_trace  <- out[[1]]   # ecs_df
fault_disp <- out[[2]]   # GAM model frame: Longitude,Latitude,wt,x,y,u,t
fit_gam_xy <- out[[3]]   # fit_nom_trace (final mgcv::gam)

# defensive: the lpmatrix anchor + GC2 truth both depend on these columns
assert_that(all(c("Longitude", "Latitude", "u", "t") %in% colnames(fault_disp)),
            msg = "fault_disp missing Longitude/Latitude/u/t; GC2 + anchor fixtures invalid.")

# --- (2) end-to-end trace + anchor + GC2 ground truth ------------------------
write.csv(ecs_trace,  file.path(out_dir, "ecs_trace.csv"),  row.names = FALSE)
write.csv(fault_disp, file.path(out_dir, "fault_disp.csv"), row.names = FALSE)

# --- (3) spline fixtures: basis (ANCHORED), coefficients, penalties ----------
lpmat <- predict(fit_gam_xy, type = "lpmatrix")
assert_that(nrow(lpmat) == nrow(fault_disp),
            msg = "lpmatrix rows != fault_disp rows; anchor broken.")
write.csv(as.data.frame(lpmat), file.path(out_dir, "lpmatrix.csv"), row.names = FALSE)

beta <- coef(fit_gam_xy)
write.csv(data.frame(name = names(beta), coef = as.numeric(beta)),
          file.path(out_dir, "coef.csv"), row.names = FALSE)

sp_used <- fit_gam_xy$sp
n_sm    <- length(fit_gam_xy$smooth)
assert_that(n_sm >= 1)
for (i in seq(n_sm)) {
  S_i <- as.matrix(fit_gam_xy$smooth[[i]]$S[[1]])
  fn  <- file(file.path(out_dir, sprintf("penalty_S_pred%d.csv", i)), "w")
  writeLines(sprintf("# sp_used=%s  lambda_p=%g  label=%s  k=%d",
                     paste(sp_used, collapse = ","), lambda_p,
                     fit_gam_xy$smooth[[i]]$label, ncol(S_i)), fn)
  write.csv(as.data.frame(S_i), fn, row.names = FALSE)
  close(fn)
}
if (n_sm == 2)
  message("penalties identical (x vs y): ",
          isTRUE(all.equal(fit_gam_xy$smooth[[1]]$S[[1]], fit_gam_xy$smooth[[2]]$S[[1]])))

message("Wrote selection, anchor/GC2, and spline fixtures to ", out_dir)
