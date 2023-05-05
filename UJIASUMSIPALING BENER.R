ara=read.csv("PAPUA.csv", sep=";", header = TRUE)
ara
x1=ara$x1
x2=ara$x2
x3=ara$x3
y=ara$y
model<-lm(y~x1+x2+x3,data = ara)
model
summary(model)
#Uji Hetero
library(lmtest)
bptest(model) #karena p-value > 0.05, maka gagal tolak H0 sehingga tidak ada hetero
#uji autokorelasi
dwtest(model) #karena pvalue > 0.05 maka terjadi autokorelasi
#uji multi
library(car)
vif(model) #karena nilai VIF < 10 maka tidak terjadi multikolinieritas
#uji normalitas
shapiro.test(residuals(model)) #karena p-value > 0.05 maka data normal