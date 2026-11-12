# Plots a SAMM run against the measurements: SOC, microbial N, microbial C, aggregate C, free MAOC, litterbag C.
# Base R only.  usage: Rscript plot_results.R output/daily_X.csv output/obs_X.csv output/X.png
args <- commandArgs(TRUE)
d <- read.csv(args[1]); m <- read.csv("AllMeasurements.csv", sep = ";", check.names = FALSE, fileEncoding = "UTF-8-BOM")
start <- as.Date("1995-03-01"); d$date <- start + d$time - 1
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
q <- list("SOC" = function(x) x$MIC_C + x$MAO_C + x$AggSTR_C + x$AggLAB_C + x$AggMAO_C,
          "Microbial N" = function(x) x$MIC_N, "Microbial C" = function(x) x$MIC_C,
          "Aggregate C" = function(x) x$AggSTR_C + x$AggLAB_C + x$AggMAO_C, "Free MOAC" = function(x) x$MAO_C,
          "Litterbag C" = function(x) x$LitterBagSTR_C + x$LitterBagLAB_C)
xlims <- list("Aggregate C" = as.Date(c("2019-01-01", "2019-12-31")), "Free MOAC" = as.Date(c("2019-01-01", "2019-12-31")),
              "Litterbag C" = as.Date(c("2007-03-01", "2008-05-31")))
png(args[3], width = 3000, height = 3400, res = 200)
par(mfrow = c(length(q), 5), mar = c(2.6, 3.6, 2.2, 0.6), mgp = c(2.2, 0.6, 0), oma = c(2, 1, 0, 0))
for (ty in names(q)) {
  xl <- if (!is.null(xlims[[ty]])) xlims[[ty]] else as.Date(c("1995-03-01", "2019-12-31"))
  mt <- m[m[["Measurement type"]] == ty, ]; mt$date <- start + mt[["Day of simulation"]] - 1
  ys <- c(q[[ty]](d)[d$date >= xl[1] & d$date <= xl[2]], mt$Value + sqrt(mt$Variance))
  for (t in names(trts)) {
    x <- d[d$TRT == t, ]; o <- mt[mt$Treatment == t, ]
    plot(x$date, q[[ty]](x), type = "l", col = "#1f77b4", lwd = 1.5, xlim = xl, ylim = c(0, max(ys, na.rm = TRUE) * 1.05),
         xlab = "", ylab = ty, main = if (ty == "SOC") trts[t] else "")
    if (nrow(o)) { arrows(o$date, pmax(o$Value - sqrt(o$Variance), 0), o$date, o$Value + sqrt(o$Variance), angle = 90,
                          code = 3, length = 0.02, col = "grey40"); points(o$date, o$Value, pch = 16, cex = 0.8) }
  }
}
mtext("Line: SAMM   Points: measured +/- 1 SD (kg/ha, 0-15 cm)", side = 1, outer = TRUE, line = 0.6)
invisible(dev.off())
