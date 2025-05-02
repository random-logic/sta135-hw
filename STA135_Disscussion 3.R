# Discussion 3

# 1. Follow the R codes to learn how to conduct multivariate normal test.

library(MVA)
library(MVN)
library(HSAUR2)
measure <- structure(list(V1 = 1:20,
                          V2 = c(34L, 37L, 38L, 36L, 38L, 43L, 40L, 38L, 40L, 41L, 36L, 36L, 34L, 33L, 36L, 37L, 34L, 36L, 38L, 35L),
                          V3 = c(30L, 32L, 30L, 33L, 29L, 32L, 33L, 30L, 30L, 32L, 24L, 25L, 24L, 22L, 26L, 26L, 25L, 26L, 28L, 23L),
                          V4 = c(32L, 37L, 36L, 39L, 33L, 38L, 42L, 40L, 37L, 39L, 35L, 37L, 37L, 34L, 38L, 37L, 38L, 37L, 40L, 35L)), .Names = c("V1", "V2", "V3", "V4"), class = "data.frame", row.names = c(NA, -20L))
measure <- measure[,-1]
names(measure) <- c("chest", "waist", "hips")
measure$gender <- gl(2, 10)
levels(measure$gender) <- c("male", "female")

toLatex(HSAURtable(measure), pcol = 2,
        caption = "Chest, waist, and hip measurements on 20 individuals (in inches).",
        label = "ch:MVA:tab:measure", rownames = FALSE)

y=matrix(0,20,3)
y[1:20,1]=measure[,"chest"];
y[1:20,2]=measure[,"waist"];
y[1:20,3]=measure[,"hips"];
qqnorm(measure[,"chest"], main = "chest"); qqline(measure[,"chest"])
qqnorm(measure[,"waist"], main = "waist"); qqline(measure[,"waist"])
qqnorm(measure[,"hips"], main = "hips"); qqline(measure[,"hips"])

ks.test(y[1:20,1],'pnorm')
ks.test(y[1:20,2],'pnorm')
ks.test(y[1:20,3],'pnorm')

shapiro.test(measure[,"chest"])
shapiro.test(measure[,"waist"])
shapiro.test(measure[,"hips"])

mvn(y,multivariatePlot = "qq")
mvn(y,mvnTest="royston")
mvn(y,mvnTest="mardia")


# 2. Testing using T^2 statistic
Y <- matrix(c(3,6,5,10,10,12,14,9),4,2)
ybar <- apply(Y,2,mean)
S <- matrix(0,2,2)
for(i in 1:4){
  S <- S + (Y[i,]-ybar)%*%t(Y[i,]-ybar)
}
S <- S/(4-1)
mu0 <- c(6,11)
T2 <- 4*(ybar-mu0) %*% solve(S) %*% (ybar-mu0)
T2

qf(0.95,2,2)
qf(0.95,2,2)*2*(4-1)/(4-2)
# not reject null hypothesis.


# 3. Real example on datasets
mydata <- read.table('nutrient.txt')
mydata <- mydata[,2:6]
Y <- mydata
ybar <- apply(Y,2,mean)
S <- cov(mydata)
n <- dim(mydata)[1]
p <- dim(mydata)[2]
#S <- matrix(0,5,5)
#for(i in 1:n){
#  S <- S + as.numeric((Y[i,]-ybar))%*%t(as.numeric(Y[i,]-ybar))
#}
#S <- S/(n-1) # the same result
mu0 <- c(1000,15,60,800,75)
T2 <- n*(ybar-mu0) %*% solve(S) %*% (ybar-mu0)
T2

qf(0.95,p,n-p)
qf(0.95,p,n-p)*p*(n-1)/(n-p)

qf(0.99,p,n-p)
qf(0.99,p,n-p)*p*(n-1)/(n-p)

# very small p-value 
# reject null hypothesis.