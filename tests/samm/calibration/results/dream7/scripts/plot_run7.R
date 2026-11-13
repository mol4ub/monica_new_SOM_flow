# Run 6: all data types; microbial C and N shown as the measured extractable quantities (no conversion of the data):
#   extractable C = MBC * k_EC (0.33, as used for the data); model: MIC_C * 0.33
#   ninhydrin-reactive N = MBN / 3.1 (as used for the data);  model: MIC_N / fN_ninhydrin (calibrated; v1.1: 3.1)
args <- commandArgs(TRUE); C <- args[1]; DATA <- args[2]; fN <- as.numeric(args[3])
kEC <- 0.33; fUsed <- 3.1
best <- read.csv(file.path(C, "dream7_daily.csv")); v11 <- read.csv(file.path(C, "v11_samm.csv"))
band <- read.csv(file.path(C, "dream7_band.csv"))
m <- read.csv(file.path(DATA, "AllMeasurements.csv"), sep = ";", check.names = FALSE, fileEncoding = "UTF-8-BOM")
start <- as.Date("1995-03-01")
for (nm in c("best", "v11", "band")) { x <- get(nm); x$date <- start + x$time - 1; assign(nm, x) }
best$ExtC <- best$MIC_C * kEC; best$NinN <- best$MIC_N / fN
v11$ExtC <- v11$MIC_C * kEC; v11$NinN <- v11$MIC_N / fUsed
# band: microbial N there is "as reported" (MIC_N * 3.1 / fN) -> ninhydrin N = / 3.1; microbial C -> * k_EC
b1 <- band$quantity == "Microbial N"; band[b1, c("q025", "q50", "q975")] <- band[b1, c("q025", "q50", "q975")] / fUsed
b2 <- band$quantity == "Microbial C"; band[b2, c("q025", "q50", "q975")] <- band[b2, c("q025", "q50", "q975")] * kEC
band$quantity[b1] <- "Ninhydrin N"; band$quantity[b2] <- "Extractable C"
# data
m$Type <- m[["Measurement type"]]; m$day <- m[["Day of simulation"]]; m$date <- start + m$day - 1; m$SD <- sqrt(m$Variance)
i <- m$Type == "Microbial N"; m$Value[i] <- m$Value[i] / fUsed; m$SD[i] <- m$SD[i] / fUsed; m$Type[i] <- "Ninhydrin N"
i <- m$Type == "Microbial C"; m$Value[i] <- m$Value[i] * kEC; m$SD[i] <- m$SD[i] * kEC; m$Type[i] <- "Extractable C"
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
q <- list("SOC" = function(x) x$LMW_C + x$MIC_C + x$MAO_C + x$AggSTR_C + x$AggLAB_C + x$AggMAO_C,
          "Ninhydrin N" = function(x) x$NinN, "Aggregate C" = function(x) x$AggSTR_C + x$AggLAB_C + x$AggMAO_C,
          "Free MOAC" = function(x) x$MAO_C, "Litterbag C" = function(x) x$LitterBagSTR_C + x$LitterBagLAB_C,
          "Extractable C" = function(x) x$ExtC,
          "Soil N" = function(x) x$LMW_N + x$MIC_N + x$MAO_N + x$AggLAB_N + x$AggMAO_N)
units <- c("kg C/ha", "kg N/ha", "kg C/ha", "kg C/ha", "kg C/ha", "kg C/ha", "kg N/ha")
xlims <- list("Aggregate C" = as.Date(c("2019-01-01", "2019-12-31")), "Free MOAC" = as.Date(c("2019-01-01", "2019-12-31")),
              "Litterbag C" = as.Date(c("2007-03-01", "2008-05-31")))
blue <- "#1f77b4"; bandcol <- adjustcolor(blue, 0.25)
drawBand <- function(b) polygon(c(b$date, rev(b$date)), c(b$q025, rev(b$q975)), col = bandcol, border = NA)
obsPts <- function(o) if (nrow(o)) { arrows(o$date, pmax(o$Value - o$SD, 0), o$date, o$Value + o$SD, angle = 90, code = 3,
                                            length = 0.02, col = "grey40"); points(o$date, o$Value, pch = 16, cex = 0.8) }
png(file.path(C, "dream7_plot.png"), width = 3000, height = 3900, res = 200)
par(mfrow = c(length(q), 5), mar = c(2.6, 3.8, 2.2, 0.6), mgp = c(2.3, 0.6, 0), oma = c(2.5, 1, 0, 0))
for (k in seq_along(q)) {
  ty <- names(q)[k]; xl <- if (!is.null(xlims[[ty]])) xlims[[ty]] else as.Date(c("1995-03-01", "2019-12-31"))
  mt <- m[m$Type == ty, ]; bt <- band[band$quantity == ty & band$date >= xl[1] & band$date <= xl[2], ]
  inr <- function(d) q[[ty]](d)[d$date >= xl[1] & d$date <= xl[2]]
  yl <- c(0, max(c(inr(best), inr(v11), bt$q975, mt$Value + mt$SD), na.rm = TRUE) * 1.05)
  for (t in names(trts)) {
    plot(NA, xlim = xl, ylim = yl, xlab = "", ylab = paste0(ty, " (", units[k], ")"), xaxt = "n", main = if (k == 1) trts[t] else "")
    if (is.null(xlims[[ty]])) axis.Date(1, at = as.Date(paste0(seq(1995, 2020, 5), "-01-01")), format = "%Y") else
      axis.Date(1, at = seq(xl[1], xl[2], by = "3 months"), format = "%b %y", cex.axis = 0.8)
    b <- bt[bt$TRT == t, ]; if (nrow(b)) drawBand(b)
    v <- v11[v11$TRT == t, ]; lines(v$date, q[[ty]](v), col = "grey60", lwd = 1)
    x <- best[best$TRT == t, ]; lines(x$date, q[[ty]](x), col = blue, lwd = 1.4)
    obsPts(mt[mt$Treatment == t, ])
  }
}
mtext(sprintf("Blue: run 7 best set (weeds, hoeing, N options), band: 95%% posterior   Grey: SAMM v1.1   Points: measured +/- 1 SD   Microbial: extractable C (= MBC x 0.33) and ninhydrin N (= MBN / 3.1); model N = MIC_N / %.2f", fN),
      side = 1, outer = TRUE, line = 1, cex = 0.9)
invisible(dev.off())

# per year with microbial data
addDays <- c(1, 382, 759, 1125, 1490, 1874, 2231, 2596, 2961, 3323, 3688, 4049, 4414, 4779, 5144, 5509, 5874, 6240,
             6605, 6970, 7335, 7702, 8068, 8432, 8797, 9162)
mm <- m[m$Type %in% c("Ninhydrin N", "Extractable C"), ]
mm$season <- sapply(mm$day, function(x) max(addDays[addDays <= x + 14])); seasons <- sort(unique(mm$season))
for (ty in c("Ninhydrin N", "Extractable C")) {
  col <- if (ty == "Ninhydrin N") "NinN" else "ExtC"
  png(file.path(C, paste0("dream7_microbes_", if (ty == "Ninhydrin N") "N" else "C", ".png")), width = 3000, height = 520 * length(seasons) + 200, res = 200)
  par(mfrow = c(length(seasons), 5), mar = c(2.4, 3.6, 1.6, 0.6), mgp = c(2.1, 0.6, 0), oma = c(2.5, 2.5, 2, 0))
  for (s in seasons) {
    obs <- mm[mm$season == s & mm$Type == ty, ]; lo <- s - 25; hi <- max(s + 250, max(mm$day[mm$season == s]) + 15)
    sel <- function(x) x[x$time >= lo & x$time <= hi, ]; bs <- sel(band[band$quantity == ty, ])
    yl <- c(0, max(c(sel(best)[[col]], sel(v11)[[col]], bs$q975, obs$Value + obs$SD), na.rm = TRUE) * 1.05)
    for (t in names(trts)) {
      x <- sel(best[best$TRT == t, ]); w <- sel(v11[v11$TRT == t, ]); o <- obs[obs$Treatment == t, ]; b <- bs[bs$TRT == t, ]
      plot(x$date, x[[col]], type = "n", ylim = yl, xlab = "", ylab = "", xaxt = "n", main = if (s == seasons[1]) trts[t] else "")
      axis.Date(1, at = seq(as.Date(format(start + lo - 1, "%Y-%m-01")), start + hi - 1, by = "2 months"), format = "%b %y", cex.axis = 0.8)
      abline(v = start + s - 1, col = "grey85", lty = 2)
      if (nrow(b)) drawBand(b)
      lines(w$date, w[[col]], col = "grey60", lwd = 1); lines(x$date, x[[col]], col = blue, lwd = 1.5); obsPts(o)
      if (t == "CT") mtext(format(start + s - 1, "%Y"), side = 2, line = 2.2, font = 2, cex = 0.9)
    }
  }
  mtext(paste0(ty, if (ty == "Ninhydrin N") " (kg N/ha; data = MBN / 3.1, model = MIC_N / fN)" else " (kg C/ha; data = MBC x 0.33, model = MIC_C x 0.33)",
               " in the years with microbial data"), side = 3, outer = TRUE, line = 0.4, cex = 1.1)
  mtext("Blue: run 7 best set, band: 95% posterior   Grey: SAMM v1.1   Points: measured +/- 1 SD   Dashed: litter addition",
        side = 1, outer = TRUE, line = 1, cex = 1)
  invisible(dev.off())
}
