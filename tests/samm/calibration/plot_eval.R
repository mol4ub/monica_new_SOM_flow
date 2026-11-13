# Evaluation plot as Fig. 5 of the SAMM paper: simulated vs measured per data type, coloured by treatment.
# Horizontal bars: measured +/- 1.96 SD; vertical bars: 95% posterior interval of the simulation (band file).
# Statistics after Loague & Green (1991) and Gauch et al. (2003): MSE = SB + NU + LC.
args <- commandArgs(TRUE); C <- args[1]; DATA <- args[2]; run <- args[3]; out <- args[4]
caption <- if (length(args) > 4) args[5] else paste0("SAMM ", run, ", best parameter set")
o <- read.csv(file.path(C, paste0(run, "_obs.csv")))
band <- read.csv(file.path(C, paste0(run, "_band.csv")))
m <- read.csv(file.path(DATA, "AllMeasurements.csv"), sep = ";", check.names = FALSE, fileEncoding = "UTF-8-BOM")
m$key <- paste(m$Treatment, m[["Day of simulation"]], m[["Measurement type"]])
o$key <- paste(o$TRT, o$day, o$type)
o$sd <- sqrt(m$Variance[match(o$key, m$key)])
bkey <- paste(band$TRT, band$time, band$quantity)
o$lo <- band$q025[match(o$key, bkey)]; o$hi <- band$q975[match(o$key, bkey)]
trts <- c(CT = "Control", DP = "Dipterocarp", GN = "Groundnut", RS = "Rice straw", TM = "Tamarind")
cols <- c(CT = "#440154", DP = "#3B528B", GN = "#21908C", RS = "#5DC863", TM = "#FDE725")
panels <- c("Aggregate C" = "Aggregate C", "Litterbag C" = "Litterbag C", "Free MOAC" = "MAO C",
            "Microbial N" = "Microbial N (as reported)", "SOC" = "SOC", "Soil N" = "Soil N", "Microbial C" = "Microbial C")
stats <- function(x, y) {
  n <- length(x); mse <- mean((y - x)^2); f <- lm(y ~ x); a <- coef(f)[1]; b <- coef(f)[2]; r2 <- summary(f)$r.squared
  sb <- (mean(y) - mean(x))^2; nu <- (1 - b)^2 * mean((x - mean(x))^2); lc <- (1 - r2) * mean((y - mean(y))^2)
  list(a = a, b = b, r2 = r2, ef = 1 - sum((y - x)^2) / sum((x - mean(x))^2), rmse = sqrt(mse),
       sb = 100 * sb / mse, nu = 100 * nu / mse, lc = 100 * lc / mse)
}
png(out, width = 3300, height = 3000, res = 220)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 8), 3, 3, byrow = TRUE))
par(mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.7, 0), cex.axis = 1.1, cex.lab = 1.2)
for (ty in names(panels)) {
  d <- o[o$type == ty & is.finite(o$sim), ]
  rng <- range(c(0, d$obs - 1.96 * d$sd, d$obs + 1.96 * d$sd, d$lo, d$hi, d$sim), na.rm = TRUE); rng[1] <- max(0, rng[1])
  plot(NA, xlim = rng, ylim = rng, xlab = "Measured value", ylab = "Simulated value", main = panels[[ty]], cex.main = 1.4)
  segments(pmax(d$obs - 1.96 * d$sd, 0), d$sim, d$obs + 1.96 * d$sd, d$sim, col = "grey75")
  segments(d$obs, d$lo, d$obs, d$hi, col = "grey75")
  abline(0, 1, lwd = 1.5)
  s <- stats(d$obs, d$sim)
  xs <- seq(min(d$obs), max(d$obs), length.out = 50); lines(xs, s$a + s$b * xs, col = "#3366FF", lwd = 2.5)
  points(d$obs, d$sim, pch = 21, bg = cols[d$TRT], col = "grey20", cex = 1.3)
  l1 <- bquote("y = " * .(sprintf("%.2f + %.2fx", s$a, s$b)) * "; " * R^2 * " = " * .(sprintf("%.2f", s$r2)))
  legend("topleft", bty = "n", cex = 0.95, legend = as.expression(c(l1,
    sprintf("EF: %.2f; RMSE: %.1f; n = %d", s$ef, s$rmse, nrow(d)),
    sprintf("SB(%%): %.0f; NU(%%): %.0f; LC(%%): %.0f", s$sb, s$nu, s$lc))))
}
plot.new(); par(mar = c(0, 0, 0, 0))
legend("top", inset = 0.15, legend = trts, pt.bg = cols, pch = 21, col = "grey20", pt.cex = 2, cex = 1.5, bty = "n", ncol = 2,
       title = "Treatment")
text(0.5, 0.35, paste0(caption, "\nblack line: 1:1, blue line: regression\n",
                       "grey bars: measured +/- 1.96 SD (horizontal),\nsimulated 95% posterior interval (vertical)\n",
                       "Microbial N as reported (model x 3.1 / fitted fN)"), cex = 1.25)
invisible(dev.off())
