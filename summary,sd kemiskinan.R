#Input Data
data=read.table("clipboard",header=T)
summary(data)
#Variabel y
var(data$y)
sd(data$y)
sd(data$x1)
sd(data$x2)
sd(data$x3)
#Menampilkan Data
View(data)

#Plot variabel x dan y
plot(data$y,data$x1)