args <- commandArgs(TRUE); C <- args[1]
d <- read.csv(file.path(C, "dream7_daily.csv")); start <- as.Date("1995-03-01"); d$date <- start + d$time - 1
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
png(file.path(C, "dream7_weeds.png"), width = 3000, height = 2000, res = 200)
par(mfrow = c(2, 1), mar = c(3, 4.5, 2, 1))
cols <- c(CT = "grey30", DP = "#8c564b", GN = "#2ca02c", RS = "#ff7f0e", TM = "#1f77b4")
for (per in list(c("2007-01-01", "2008-12-31"), c("2018-01-01", "2019-12-31"))) {
  xl <- as.Date(per); s <- d[d$date >= xl[1] & d$date <= xl[2], ]
  plot(NA, xlim = xl, ylim = c(0, max(s$WeedAGB / 1000, 2.3)), xlab = "", ylab = "weed above-ground biomass (t DM/ha)", xaxt = "n",
       main = paste("Simulated weed biomass", format(xl[1], "%Y"), "-", format(xl[2], "%Y")))
  axis.Date(1, at = seq(xl[1], xl[2], by = "2 months"), format = "%b %y", cex.axis = 0.8)
  rect(as.Date(paste0(format(xl[1], "%Y"), "-04-01")) + c(0, 365), 0, as.Date(paste0(format(xl[1], "%Y"), "-09-30")) + c(0, 365), 10, col = "grey95", border = NA)
  for (t in names(trts)) { x <- s[s$TRT == t, ]; lines(x$date, x$WeedAGB / 1000, col = cols[t], lwd = 1.5) }
  abline(h = c(1.5, 0.3), lty = 2, col = c(cols["GN"], cols["CT"]))
  legend("topleft", c(trts, "your estimate GN (1.5) / CT (0.3)"), col = c(cols, "grey50"), lty = c(rep(1, 5), 2), lwd = 1.5, bty = "n", cex = 0.8)
}
invisible(dev.off())
