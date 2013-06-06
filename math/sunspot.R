setwd("~/research/var/sunspot/expr/001-initial")
R = rbind(load.files(find.files("ta0", "fitness.dat.gz")))
R$normalized_fitness = 100.0 / R$max_fitness - 1.0

quartz(width=6,height=3.75)
g = ggplot(data=subset(R,update%%100==0), aes(x=update, y=normalized_fitness)) 
g + geom_line(aes(color=trial))


setwd("~/research/var/sunspot/expr/001-initial")
setwd("/Users/dk/research/src/sunspot")

D = rbind(load.files(find.files("tmp", "test_sunspot.dat"),tag=FALSE))

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red")
plot(D$predicted, type="l",col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red", xlim=c(3900,4100))
lines(D$predicted, col="blue")

quartz(width=6,height=3.75)
plot(D$observed, type="l",col="red", xlim=c(6000,6200))
lines(D$predicted, col="blue")

D$diff = D$observed-D$predicted
quartz(width=6,height=3.75)
g = ggplot(data=D,aes(x=1:nrow(D),y=diff))
g + geom_line()

quartz(width=6,height=3.75)
g = ggplot(data=D,aes(x=observed, y=predicted))
g + geom_point()