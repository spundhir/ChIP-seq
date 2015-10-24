## help(distribution) to see other distributions
ARGV <- commandArgs(TRUE)
pnorm(as.numeric(ARGV[1]), as.numeric(ARGV[2]), as.numeric(ARGV[3]), lower.tail=F)
q()
