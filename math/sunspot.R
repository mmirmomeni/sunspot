
### fitness
setwd("~/research/var/sunspot/002-multi")
R = rbind(load.files(find.files("ta0", "fitness.dat.gz")))
R$normalized_fitness = 100.0 / R$max_fitness - 1.0

quartz(width=6,height=3.75)
g = ggplot(data=subset(R,update%%100==0), aes(x=update, y=normalized_fitness)) 
g + geom_line(aes(color=trial))


### test results
setwd("/Users/dk/research/src/sunspot")
D = rbind(load.files(find.files("tmp", "test_sunspot.dat"),tag=FALSE))
l = length(D$observed) - 1

### deltas
quartz()
par(mfrow=c(4,1))
offset=1
d1 = D$predicted_tplus1[0:-offset] - D$observed[offset:l]
plot(d1, type="l", col="red")
title("(t+1) - observed")

offset=2
d2 = D$predicted_tplus2[0:-offset] - D$observed[offset:l]
#quartz(width=6,height=3.75)
plot(d2, type="l", col="red")
title("(t+2) - observed")

offset=3
d3 = D$predicted_tplus3[0:-offset] - D$observed[offset:l]
#quartz(width=6,height=3.75)
plot(d3, type="l", col="red")
title("(t+3) - observed")

offset=4
d4 = D$predicted_tplus4[0:-offset] - D$observed[offset:l]
#quartz(width=6,height=3.75)
plot(d4, type="l", col="red")
title("(t+4) - observed")

### individual zoom
quartz(width=6,height=3.75)
offset=1
plot(D$observed[offset:l], type="l",col="blue", xlim=c(3900,4100))
lines(D$predicted_tplus1[0:-offset], col="red")
legend("topleft", c("observed","t+1"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))

quartz(width=6,height=3.75)
offset=2
plot(D$observed[offset:l], type="l",col="blue", xlim=c(3900,4100))
lines(D$predicted_tplus2[0:-offset], col="red")
legend("topleft", c("observed","t+2"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))

quartz(width=6,height=3.75)
offset=3
plot(D$observed[offset:l], type="l",col="blue", xlim=c(3900,4100))
lines(D$predicted_tplus3[0:-offset], col="red")
legend("topleft", c("observed","t+3"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))

quartz(width=6,height=3.75)
offset=4
plot(D$observed[offset:l], type="l",col="blue", xlim=c(3900,4100))
lines(D$predicted_tplus4[0:-offset], col="red")
legend("topleft", c("observed","t+4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"))


quartz(width=6,height=3.75)
plot(D$observed[1:l], type="l",col="red")
plot(D$predicted_tplus1[0:-1], type="l",col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red")
plot(D$predicted_tplus2, type="l",col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red")
plot(D$predicted_tplus3, type="l",col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red")
plot(D$predicted_tplus4, type="l",col="blue")



quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red", xlim=c(3900,4100))
lines(D$predicted_tplus4, col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red", xlim=c(6000,6200))
lines(D$predicted, col="blue")

D$diff = D$observed-D$predicted_tplus4
quartz(width=6,height=3.75)
g = ggplot(data=D,aes(x=1:nrow(D),y=diff))
g + geom_line()

quartz(width=6,height=3.75)
g = ggplot(data=D,aes(x=observed, y=predicted))
g + geom_point()