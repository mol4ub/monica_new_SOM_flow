# Plots the MONICA-SAMM results of the Khon Kaen experiment against the measurements and prints RMSE and EF.
# Base R only. Run in this folder after run_khonkaen.bat:  Rscript plot_results.R
# Writes output/khonkaen_results.png and output/matched.csv.
start <- as.Date("1995-03-01")   # day 1 of the simulation (measurements are given by "Day of simulation")
kgha <- 0.15 * 10000             # kg m-3 in the 0-15 cm layer -> kg ha-1
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
types <- c("SOC", "Soil N", "Microbial C", "Microbial N", "Aggregate C", "Free MOAC")
units <- c("kg C/ha", "kg N/ha", "kg C/ha", "kg N/ha", "kg C/ha", "kg C/ha")

read_monica <- function(f) {
  hdr <- strsplit(readLines(f, n = 2)[2], ",")[[1]]
  d <- read.csv(f, skip = 3, header = FALSE, col.names = hdr)
  p <- function(n) d[[paste0("SAMM_", n, "_1")]] * kgha
  data.frame(date = as.Date(d$Date),
             "SOC" = p("MIC_C") + p("MAO_C") + p("AggSTR_C") + p("AggLAB_C") + p("AggMAO_C"),
             "Soil N" = p("LMW_N") + p("MIC_N") + p("MAO_N") + p("AggLAB_N") + p("AggMAO_N"),
             "Microbial C" = p("MIC_C"), "Microbial N" = p("MIC_N"),
             "Aggregate C" = p("AggSTR_C") + p("AggLAB_C") + p("AggMAO_C"), "Free MOAC" = p("MAO_C"),
             check.names = FALSE)
}
sim <- list()
for (t in names(trts)) {
  f <- file.path("output", paste0("out-", t, ".csv"))
  if (file.exists(f)) sim[[t]] <- read_monica(f)
}
if (length(sim) == 0) stop("no output/out-<treatment>.csv found; run run_khonkaen.bat first")

m <- read.csv("measurements.csv", sep = ";", check.names = FALSE, fileEncoding = "UTF-8-BOM")
m <- m[m[["Measurement type"]] %in% types & m$Treatment %in% names(sim), ]
# MONICA writes the state at the end of a day; measurement day t = state at the start of day t
m$date <- start + pmax(m[["Day of simulation"]], 2) - 2
m$sd <- sqrt(m$Variance)
m$sim <- mapply(function(t, ty, d) { s <- sim[[t]]; v <- s[[ty]][s$date == d]; if (length(v)) v else NA },
                m$Treatment, m[["Measurement type"]], m$date)
write.csv(m, file.path("output", "matched.csv"), row.names = FALSE)

cat(sprintf("%-12s %5s %9s %7s %9s\n", "type", "n", "RMSE", "EF", "bias"))
for (ty in types) {
  x <- m[m[["Measurement type"]] == ty & !is.na(m$sim), ]
  if (nrow(x) == 0) next
  rmse <- sqrt(mean((x$sim - x$Value)^2))
  ef <- 1 - sum((x$sim - x$Value)^2) / sum((x$Value - mean(x$Value))^2)
  cat(sprintf("%-12s %5d %9.1f %7.2f %+9.1f\n", ty, nrow(x), rmse, ef, mean(x$sim - x$Value)))
}

png(file.path("output", "khonkaen_results.png"), width = 3000, height = 3400, res = 200)
par(mfrow = c(length(types), length(sim)), mar = c(2.6, 3.6, 2.2, 0.6), mgp = c(2.2, 0.6, 0), oma = c(2, 1, 0, 0))
for (i in seq_along(types)) {
  ty <- types[i]
  zoom <- ty %in% c("Aggregate C", "Free MOAC")   # measured only in 2019
  xl <- if (zoom) as.Date(c("2019-01-01", "2019-12-31")) else as.Date(c("1995-03-01", "2019-12-31"))
  mt <- m[m[["Measurement type"]] == ty, ]
  yl <- c(0, 1.05 * max(c(unlist(lapply(sim, function(s) s[[ty]][s$date >= xl[1] & s$date <= xl[2]])),
                          mt$Value + mt$sd), na.rm = TRUE))
  for (t in names(sim)) {
    s <- sim[[t]]; x <- mt[mt$Treatment == t, ]
    plot(s$date, s[[ty]], type = "l", col = "#1f77b4", lwd = 1.5, xlim = xl, ylim = yl, xlab = "",
         ylab = paste0(ty, " (", units[i], ")"), main = if (i == 1) trts[t] else "")
    if (nrow(x) > 0) {
      arrows(x$date, pmax(x$Value - x$sd, 0), x$date, x$Value + x$sd, angle = 90, code = 3, length = 0.02, col = "grey30")
      points(x$date, x$Value, pch = 16, cex = 0.9)
    }
  }
}
mtext("Line: SAMM in MONICA    Points: measured +/- 1 SD", side = 1, outer = TRUE, line = 0.6, cex = 1.1)
invisible(dev.off())
cat("written: output/khonkaen_results.png, output/matched.csv\n")
