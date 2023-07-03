## nth order FDN reverbrater
library(tuneR)
library(signal)
library(numbers)
cwd <- getwd()

source(paste(cwd, "had_matrix.R", sep = "/"))
# set coefficients
orders <- 8
delaytime <- c(13, 49, 71, 137, 141, 149, 137, 141)
feedback <- 0.7
cutoff <- 1000
B <- matrix(rep(1, orders), nrow = orders)
C <- B * 0.7
R <- hadamard(orders)
G <- R * 0.58

# set file
prefix <- paste(cwd, "input", sep = "/")
fnames <- list.files(path = prefix, pattern = "*.wav")

for (fname in fnames) {
  setwd(paste(cwd, "input", sep = "/"))
  wav <- fname
  input <- readWave(wav)@left
  srate <- readWave(wav)@samp.rate

  # coefficients for LPF
  fn <- srate / 2
  ff <- cutoff / fn
  fp <- ff * pi
  b1 <- 1
  b2 <- 0.12
  a1 <- 1
  a2 <- exp((-1) * fp)
  norm <- (1 - a2) / (b1 + b2)
  b1 <- b1 * norm
  b2 <- b2 * norm

  # LPF
  signal <- numeric(length(input))
  n <- 1
  while (n <= length(input)) {
    if (n == 1) {
      signal[n] <- input[n]
      n <- n + 1
    } else {
      outb1 <- input[n] * b1
      outb2 <- input[(n - 1)] * b2
      signal[n] <- outb1 + outb2 + (signal[(n - 1)] * a2)
      n <- n + 1
    }
  }
  n <- 1
  # ring buffer
  buffer <- matrix(numeric(ceiling(max(srate * delaytime / 1000) * orders)), ncol = orders)
  bf_w <- matrix(round(srate * delaytime / 1000), ncol = orders)
  bf_r <- matrix(rep(1, orders), nrow = orders)

  # calculation
  out <- numeric(length(signal))
  ## read
  nrb <- 1
  for (n in 1:length(signal)) {
    bf_out <- matrix(numeric(orders), ncol = 1)
    for (nrb in 1:orders) {
      bf_out[nrb, ] <- buffer[bf_r[nrb, ], nrb]
    }
    out[n] <- sum(bf_out * C)
    ## write
    bf_in <- signal[n] * B
    for (nrb in 1:orders) {
      buffer[bf_w[, nrb], nrb] <- bf_in[nrb, ] + (sum(bf_out * G[, nrb])) * feedback
    }
    for (nrb in 1:orders) {
      bf_r[nrb, ] <- bf_r[nrb, ] %% nrow(buffer) + 1
      bf_w[, nrb] <- bf_w[, nrb] %% nrow(buffer) + 1
    }
  }
  # save wav
  dat <- Wave(left = out, right = out, samp.rate = srate, bit = 32, pcm = TRUE)
  setwd(paste(cwd, "output/matrix", sep = "/"))
  writeWave(normalize(dat, unit = "32", center = TRUE), filename = sprintf("%s_%02dth_%05d.wav", substring(fname, 1, (nchar(fname) - 4)), orders, cutoff))
}
