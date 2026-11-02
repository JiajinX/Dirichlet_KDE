invisible(lapply(dir("R", full.names = TRUE), source))

reps <- 100 # Number of repetitions
r <- 0 # 0,1: density / first derivative
n <- 100 # 100,300,500,700: sample size
h <- 0.01 # bandwidth: 0.01 for r=0; 0.015 for r=1
mix0 <- Dir_mix(2)

# Exact MISE representation
S <- 5e5
mise <- time <- numeric(reps)
for (i in seq_len(reps)) {
  time[i] <- system.time({
    mix0 <- gene_compsamp(mix0, S)
    mise[i] <- if (r == 0) MISEf(h, mix0, n) else DMISEf(h, mix0, n)
  })["elapsed"]
}

# Direct Monte Carlo approximation
B <- 1e2
Q <- 1e3
mise <- time <- numeric(reps)
for (i in seq_len(reps)) {
  time[i] <- system.time({
    mix0 <- gene_data(n, mix0, B)
    mix0 <- gene_eval(mix0, neva = Q, r = r)
    mise[i] <- MISEmc(h, mix0, r = r)
  })["elapsed"]
}

# Results are saved in: results/MISE_time_SE/
