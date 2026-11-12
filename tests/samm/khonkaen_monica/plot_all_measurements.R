# Simulated (SAMM in MONICA, default MONICA modifiers; SAMM v1.1 in R as reference) vs measured, for every
# measurement type and treatment of the Khon Kaen experiment.
# Input: all_series.csv and all_obs.csv of all_measurements.py.
# usage: Rscript plot_all_measurements.R <dir> <out.png>
args <- commandArgs(TRUE)
s <- read.csv(file.path(args[1], "all_series.csv"))
o <- read.csv(file.path(args[1], "all_obs.csv"))
s$date <- as.Date(s$date); o$date <- as.Date(o$date)
types <- c("SOC", "Soil N", "Microbial C", "Microbial N", "Aggregate C", "Free MOAC")
units <- c("kg C/ha", "kg N/ha", "kg C/ha", "kg N/ha", "kg C/ha", "kg C/ha")
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
xlims <- list("Aggregate C" = as.Date(c("2019-01-01", "2019-12-31")),
              "Free MOAC" = as.Date(c("2019-01-01", "2019-12-31")))
blue <- "#1f77b4"; grey <- "grey55"
png(args[2], width = 3000, height = 3400, res = 200)
par(mfrow = c(length(types), length(trts)), mar = c(2.6, 3.6, 2.2, 0.6), mgp = c(2.2, 0.6, 0), oma = c(3, 1, 0, 0))
for (i in seq_along(types)) {
  ty <- types[i]
  xl <- if (!is.null(xlims[[ty]])) xlims[[ty]] else as.Date(c("1995-03-01", "2019-12-31"))
  st <- s[s$type == ty & s$date >= xl[1] & s$date <= xl[2], ]
  ot <- o[o$type == ty, ]
  vals <- c(suppressWarnings(as.numeric(st$MONICA)), suppressWarnings(as.numeric(st$R)), ot$obs + ot$sd)
  yl <- c(0, max(vals, na.rm = TRUE) * 1.05)
  for (trt in names(trts)) {
    x <- st[st$TRT == trt, ]
    m <- ot[ot$TRT == trt, ]
    plot(NA, xlim = xl, ylim = yl, xlab = "", ylab = paste0(ty, " (", units[i], ")"), xaxt = "n",
         main = if (i == 1) trts[trt] else "", cex.main = 1.3)
    if (ty %in% c("Aggregate C", "Free MOAC")) {
      axis.Date(1, at = seq(xl[1], xl[2], by = "2 months"), format = "%b %y", cex.axis = 0.8)
    } else {
      axis.Date(1, at = as.Date(paste0(seq(1995, 2020, 5), "-01-01")), format = "%Y")
    }
    if (nrow(x) == 0 && nrow(m) == 0) { text(mean(xl), mean(yl), "no data", col = grey); next }
    lines(x$date, suppressWarnings(as.numeric(x$R)), col = grey, lwd = 1)
    lines(x$date, suppressWarnings(as.numeric(x$MONICA)), col = blue, lwd = 1.6)
    if (nrow(m) > 0) {
      arrows(m$date, pmax(m$obs - m$sd, 0), m$date, m$obs + m$sd, angle = 90, code = 3, length = 0.02, col = "grey30")
      points(m$date, m$obs, pch = 16, cex = 0.9, col = "black")
    }
  }
}
mtext("Blue: SAMM in MONICA (MONICA modifiers, default)    Grey: SAMM v1.1 in R (HYDRUS drivers)    Points: measured +/- 1 SD",
      side = 1, outer = TRUE, line = 1.2, cex = 1.1)
dev.off()
