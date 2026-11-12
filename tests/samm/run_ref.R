# Runs SAMM v1.1 (tests/samm/data/Samm_main.R, without its plotting part) and writes all treatments to a CSV.
# usage (in tests/samm/data, needs the R package deSolve): Rscript ../run_ref.R Samm_main.R R_out.csv
args <- commandArgs(TRUE); src <- args[1]; out <- args[2]
lines <- readLines(src)
end <- grep("^outComb\\$AgC", lines) - 1          # stop before the plotting part
lines <- lines[1:end]
lines <- lines[!grepl("^library\\(\"ggplot2\"\\)|^rm\\(|^cat\\(", lines)]
eval(parse(text = lines))
write.csv(outComb, out, row.names = FALSE)
cat("rows", nrow(outComb), "\n")
