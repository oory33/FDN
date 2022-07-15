library(tuneR)
library(signal)
library(numbers)
setwd("/Users/ryo/Documents/R/R_study/FDN/input")


#set file
prefix <- "/Users/ryo/Documents/R/R_study/FDN/input/"
fnames <- list.files(path=prefix, pattern="*.wav")


for (fname in fnames) {
wav <- fname
signal <- readWave(wav)@left
srate <- readWave(wav)@samp.rate
delay <- 200
feedback <- 0.6
out <- numeric(length(signal))

##buffer length
bufflen <- srate * delay / 1000

#output
n <- 1
while(n <= length(signal)){
  if(n <= bufflen){
    out[n] <- signal[n]
    n = n + 1
  } else {
    out[n:n] <- signal[n] +
    (out[(n-bufflen)] * feedback)
    n = n + 1
  }
}

dat <- Wave(left=out, right=out, samp.rate=48000, bit=32, pcm=TRUE)
setwd("/Users/ryo/Documents/R/R_study/output/fd_out")
writeWave(normalize(dat, unit="32", center=TRUE), filename=sprintf("delay_%s_%03d.wav",  substring(fname, 1, (nchar(fname)- 4)), delay))
setwd("/Users/ryo/Documents/R/R_study/FDN/input")
}
