##8th order FDN reverbrater
library(tuneR)
library(signal)
library(numbers)
cwd <- getwd()

#set file
prefix <- paste(cwd, "input", sep = "/")
fnames <- list.files(path=prefix, pattern="*.wav")

for(fname in fnames){
 setwd(paste(cwd, "input", sep = "/"))
wav <- fname
signal <- readWave(wav)@left
srate <- readWave(wav)@samp.rate
##set parameter
t066d <- c(5, 13, 37, 47, 53, 71, 83, 89)
t066g <- c(0.8, 0.6, 0.6, 0.5, 0.5, 0.5, 0.3, 0.3)
qs066 <- 0.54
t107d <- c(5, 13, 73, 97, 101, 107, 157, 151)
t107g <- c(0.8, 0.6, 0.6, 0.5, 0.5, 0.5, 0.3, 0.3)
qs107 <- 0.5
t23d <- c(13, 49, 71, 137, 141, 149, 137, 149)
t23g <- c(0.8, 0.6, 0.6, 0.5, 0.5, 0.5, 0.3, 0.3)
qs23 <- 0.57
t423d <- c(5, 13, 73, 97, 151, 157, 307, 311)
t423g      <- c(0.8, 0.6, 0.6, 0.5, 0.5, 0.5, 0.3, 0.3)
qs423 <- 0.65

delaytimes <-t423d
gains      <- t423g
qs <- 0.688
c1_delay <- delaytimes[1]
c2_delay <- delaytimes[2]
c3_delay <- delaytimes[3]
c4_delay <- delaytimes[4]
c5_delay <- delaytimes[5]
c6_delay <- delaytimes[6]
c7_delay <- delaytimes[7]
c8_delay <- delaytimes[8]
c1_gain <- gains[1]
c2_gain <- gains[2]
c3_gain <- gains[3]
c4_gain <- gains[4]
c5_gain <- gains[5]
c6_gain <- gains[6]
c7_gain <- gains[7]
c8_gain <- gains[8]
cutoff <- 1000
#pp for direct sound
mix <- 0.64

##convert coefficients
d1 <- srate * c1_delay / 1000
d2 <- srate * c2_delay / 1000
d3 <- srate * c3_delay / 1000
d4 <- srate * c4_delay / 1000
d5 <- srate * c5_delay / 1000
d6 <- srate * c6_delay / 1000
d7 <- srate * c7_delay / 1000
d8 <- srate * c8_delay / 1000
g1 <- c1_gain
g2 <- c2_gain
g3 <- c3_gain
g4 <- c4_gain
g5 <- c5_gain
g6 <- c6_gain
g7 <- c7_gain
g8 <- c8_gain

##buffer for maimum delay time
bufflength <- max(d1, d2, d3, d4, d5, d6, d7, d8)
min <- min(d1, d2, d3, d4, d5, d6, d7, d8)

#make buffer
out <- numeric(length(signal))
buf1 <- c(out,out)
buf2 <- buf1
buf3 <- buf1
buf4 <- buf1
buf5 <- buf1
buf6 <- buf1
buf7 <- buf1
buf8 <- buf1

##coefficients of LPF
fn = srate / 2
ff = cutoff / fn
fp = ff * pi
b1 = 1
b2 = 0.12
a1 = 1
a2 = exp((-1) * fp)
norm = (1 - a2) / (b1 + b2)
b1 = b1 * norm
b2 = b2 * norm

##coefficients for matrix
gp = 1 / sqrt(2)
gn = -1 * gp
q = qs

##LPF
n <- 1
while(n <= length(signal)) {
  if(n == 1){
    out[n] <- signal[n]
    n = n + 1
  } else {
 outb1 <- signal[n] * b1
 outb2 <- signal[(n-1)] * b2
 out[n] <- outb1 + outb2 + (out[(n-1)] * a2)
 n = n + 1
  }
}

#fit zero to head
sig <- c(numeric(bufflength), out)
preout <- c(numeric(length(sig)))

##FDN
n <- bufflength + 1
while(n <= length(sig) + min){
  preout[n] <- (((buf1[(n-d1)] * g1) + (buf2[(n-d2)] * g2) +
                 (buf3[(n-d3)] * g3) + (buf4[(n-d4)] * g4) +
                 (buf5[(n-d5)] * g5) + (buf8[(n-d6)] * g6) +
                 (buf7[(n-d7)] * g7) + (buf8[(n-d8)] * g8)) / 8) * (1 - mix) + sig[(n - min)] * mix
  buf1[n] <- sig[n] + ((buf2[(n-d2)] * gp) + (buf3[(n-d3)] * gp) + (buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gp) +
                       (buf5[(n-d5)] * gp) + (buf6[(n-d6)] * gp) + (buf7[(n-d7)] * gp) + (buf8[(n-d8)] * gp)) * q
  buf2[n] <- sig[n] + ((buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gn) + (buf2[(n-d2)] * gn) + (buf3[(n-d3)] * gp) +
                       (buf5[(n-d5)] * gp) + (buf6[(n-d6)] * gn) + (buf7[(n-d7)] * gp) + (buf8[(n-d8)] * gn)) * q
  buf3[n] <- sig[n] + ((buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gn) + (buf2[(n-d2)] * gp) + (buf3[(n-d3)] * gn) +
                       (buf5[(n-d5)] * gp) + (buf6[(n-d6)] * gp) + (buf7[(n-d7)] * gn) + (buf8[(n-d8)] * gn)) * q
  buf4[n] <- sig[n] + ((buf2[(n-d2)] * gn) + (buf3[(n-d3)] * gn) + (buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gp) +
                       (buf5[(n-d5)] * gp) + (buf6[(n-d6)] * gn) + (buf7[(n-d7)] * gn) + (buf8[(n-d8)] * gp)) * q
  buf1[n] <- sig[n] + ((buf2[(n-d2)] * gp) + (buf3[(n-d3)] * gp) + (buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gp) +
                       (buf5[(n-d5)] * gn) + (buf6[(n-d6)] * gn) + (buf7[(n-d7)] * gn) + (buf8[(n-d8)] * gn)) * q
  buf2[n] <- sig[n] + ((buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gn) + (buf2[(n-d2)] * gn) + (buf3[(n-d3)] * gp) +
                       (buf5[(n-d5)] * gn) + (buf6[(n-d6)] * gp) + (buf7[(n-d7)] * gn) + (buf8[(n-d8)] * gp)) * q
  buf3[n] <- sig[n] + ((buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gn) + (buf2[(n-d2)] * gp) + (buf3[(n-d3)] * gn) +
                       (buf5[(n-d5)] * gn) + (buf6[(n-d6)] * gn) + (buf7[(n-d7)] * gp) + (buf8[(n-d8)] * gp)) * q
  buf4[n] <- sig[n] + ((buf2[(n-d2)] * gn) + (buf3[(n-d3)] * gn) + (buf1[(n-d1)] * gp) + (buf4[(n-d4)] * gp) +
                       (buf5[(n-d5)] * gn) + (buf6[(n-d6)] * gp) + (buf7[(n-d7)] * gp) + (buf8[(n-d8)] * gn)) * q
  n = n + 1
}
output <- preout[(bufflength + min + 1):(length(sig) + min)]

#save wav
 dat <- Wave(left=output, right=output, samp.rate=48000, bit=32, pcm=TRUE)
 setwd(paste(cwd, "output/matrix", sep = "/"))
 writeWave(normalize(dat, unit="32", center=TRUE), filename=sprintf("%s_%05d.wav",  substring(fname, 1, (nchar(fname)- 4)), cutoff))
}
