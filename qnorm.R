## help(distribution) to see other distributions
ARGV <- commandArgs(TRUE)
qnorm(0.05, as.numeric(ARGV[1]), as.numeric(ARGV[2]),lower.tail = F, log.p = FALSE)
q()
