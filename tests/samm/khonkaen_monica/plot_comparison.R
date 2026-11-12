# Plots SOC (MIC + MAO + aggregate pools, 0-15 cm) of SAMM v1.1 in R and SAMM in MONICA for the Khon Kaen
# treatments, with the measurements. Input: comparison.csv and measured_soc.csv of compare_monica_r.py.
# usage: Rscript plot_comparison.R <dir with comparison.csv> <out.png>
args <- commandArgs(TRUE)
d <- read.csv(file.path(args[1], "comparison.csv"))
m <- read.csv(file.path(args[1], "measured_soc.csv"))
names <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
cols <- c(R = "black", MONICA = "#1f77b4", SAMM = "#d62728")
png(args[2], width = 2400, height = 1500, res = 200)
par(mfrow = c(2, 3), mar = c(4, 4.5, 2.5, 1))
ylim <- c(0, max(c(d$SOC_R, d$SOC_SAMM, d$SOC_MONICA, m$SOC), na.rm = TRUE) * 1.02)
for (trt in names(names)) {
  x <- d[d$TRT == trt, ]
  yr <- 1995 + (x$time - 1) / 365.25
  plot(yr, x$SOC_R, type = "l", lwd = 2, col = cols["R"], ylim = ylim, xlab = "Year",
       ylab = expression(SOC ~ (kg ~ C ~ ha^{-1} ~ "0-15 cm")), main = names[trt])
  lines(yr, x$SOC_MONICA, lwd = 2, col = cols["MONICA"])
  lines(yr, x$SOC_SAMM, lwd = 1.5, col = cols["SAMM"], lty = 2)
  mm <- m[m$TRT == trt, ]
  points(1995 + (mm$time - 1) / 365.25, mm$SOC, pch = 16, cex = 0.9, col = "grey35")
}
plot.new()
legend("center", bty = "n", cex = 1.05, lwd = c(2, 2, 1.5, NA), lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16),
       col = c(cols["R"], cols["MONICA"], cols["SAMM"], "grey35"),
       legend = c("SAMM v1.1 in R\n(HYDRUS soil drivers)",
                  "MONICA + SAMM,\nMONICA modifiers (default)",
                  "MONICA + SAMM, SAMM modifiers\n(pF from MONICA's retention curve)",
                  "Measured"), y.intersp = 1.6)
dev.off()
