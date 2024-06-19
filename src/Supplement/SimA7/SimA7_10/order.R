library(R.matlab)
library("lbfgs")
library("survPresmooth")
library("MASS")
library("mvtnorm")
library("survival")
library("ncvreg")
library("doParallel")
library("foreach")
library("glmnet")
source("Gen_Dat.R")
source("Fun.R")
source("screen.R")
source("method.R")
# 设置为该函数所在的路径
setwd("C:/Users/dell/Desktop/2023.0257/scr/Supplement/SimA7/SimuA7_10")

# 输入参数
inipara <- readMat("para.mat");


X = inipara$x.sc
Y = inipara$y.sc

Delta = inipara$Delta.sc

stat.scr = IPOD(X,Delta,Y,gamma=0.8)
or = c(order(stat.scr, decreasing = T))


writeMat("screen.mat",r_ans=or)
