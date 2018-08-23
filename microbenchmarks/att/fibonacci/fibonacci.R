
a <- 0
b <- 0
phi <- 1.6180339887498949

a <- floor(runif(3500000)*1000)
b <- (phi^a - (-phi)^(-a))/sqrt(5)

# We need to do something with the result
# otherwise Renjin will optimize it away

print(mean(b))
