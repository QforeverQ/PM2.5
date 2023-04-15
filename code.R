library(ggplot2)
library(gridExtra)
library(glmnet)
library(car)   # 模型诊断
library(mgcv)  # 广义可加模型

setwd("D:/研究生/统计基础/期末大作业")
# 导入训练集和测试集数据
data_train <- read.csv("D:/研究生/统计基础/期末大作业/西安pm2.5数据训练集.csv",header=TRUE,fileEncoding='UTF-8')
data_test <- read.csv("D:/研究生/统计基础/期末大作业/西安pm2.5数据测试集.csv",header=TRUE)
#将数据变为数据框
data_train<-as.data.frame(data_train)
data_test<-as.data.frame(data_test)

head(data_train)
attach(data_train)
# 查看数据结构
str(data_train)
# 查看数据的概览信息
summary(data_train)

# PM2.5与15个协变量间的散点图
for (i in 1:15) {
  pg <- ggplot(data = data_train, mapping = aes(x = data_train[,i+2], y = data_train[,2])) + 
    geom_point(fill = 'steelblue') + labs(x = names(data_train)[i+2], y = 'PM2.5')
  ggsave(paste0(i, "_", ".png"), pg, width = 3, height = 3)
}

# 建模：多元线性回归模型
liner = lm(PM2.5 ~ PM10 + SO2 + NO2 + CO + O3 + MAXT + AT + MINT + AH + MAXP + AP +
             MINP + MAXW + RAIN + MAXR, data = data_train)
summary(liner)

# 逐步回归完成遍历选择
liner_step <- step(liner)
summary(liner_step)


# 模型诊断
# 正态性检验
# 图形法
#绘制直方图
hist(x = data_train$PM2.5, freq = FALSE, main = 'PM2.5的直方图',
     ylab = '核密度值',xlab = NULL, ylim = c(0,0.02), col = 'steelblue')
#添加核密度图
lines(density(data_train$PM2.5), col = 'red', lty = 1, lwd = 2)
#添加正态分布图
x <- data_train$PM2.5[order(data_train$PM2.5)]
lines(x, dnorm(x, mean(x), sd(x)),
      col = 'black', lty = 2, lwd = 2.5)

#添加图例
legend('topright',legend = c('核密度曲线','正态分布曲线'),
       col = c('red','black'), lty = c(1,2),
       lwd = c(2,2.5), bty = 'n')

# PP图
real_dist <- ppoints(data_train$PM2.5)
theory_dist <- pnorm(data_train$PM2.5, mean = mean(data_train$PM2.5), 
                     sd = sd(data_train$PM2.5))

# 绘图
plot(sort(theory_dist), real_dist, col = 'steelblue', 
     pch = 20, main = 'PP图', xlab = '理论累计概率', 
     ylab = '实际累计概率')

# 添加对角线作为参考线
abline(a = 0,b = 1, col = 'red', lwd = 2)

# QQ图
qqnorm(data_train$PM2.5, col = 'steelblue', pch = 20,
       main = 'QQ图', xlab = '理论分位数', 
       ylab = '实际分位数')
# 绘制参考线
qqline(data_train$PM2.5, col = 'red', lwd = 2)

# 统计法
shapiro.test(data_train$PM2.5)

# 当PM2.5不满足正态分布时
# box-cox变换
powerTransform(liner_step)   
#开根号
liner_sqrt <- lm(sqrt(PM2.5) ~ PM10 + SO2 + NO2 + CO + O3 + MAXT + AT + MINT + AH + MAXP + AP +
                  MINP + MAXW + RAIN + MAXR, data = data_train)
summary(liner_sqrt)
liner_sqrt_step <- step(liner_sqrt)
summary(liner_sqrt_step)

# 多重共线性检验
vif(liner_sqrt_step)  # 通过
plot(liner_sqrt_step)



# 方差齐性检验：不通过
ncvTest(liner_sqrt_step)

# 加权最小二乘法
W<-diag(1/abs(liner_sqrt_step$residuals))
M<-as.matrix(data_train)
M = M[,2:17]
M=apply(M,2,as.numeric)
WM <- W %*% M
data_train_W <- as.data.frame(WM)

liner_weight <- lm(PM2.5 ~ PM10 + SO2 + NO2 + CO + O3 + MAXT + AT + MINT + AH + MAXP + AP +
                         MINP + MAXW + RAIN + MAXR, data = data_train_W)
summary(liner_weight)
liner_weight_step <- step(liner_weight)
summary(liner_weight_step)
powerTransform(liner_weight_step)  #不需要BOX-COX变换
vif(liner_weight_step)  #有共曲线性

# 根据变量显著性和共曲线性进行变量筛选
liner_weight_vif <- lm(PM2.5 ~ PM10 + NO2 + AH + MAXR, data = data_train_W)
summary(liner_weight_vif)
vif(liner_weight_vif)

plot(liner_weight_vif)



#数据转化为矩阵
M = as.matrix(data_train)
X = M[,3:17]
X=apply(X,2,as.numeric)
y = M[,2]

#岭回归：效果好
#可视化确定lambda
ridge<-glmnet(X,y,alpha=0)
coef(ridge)
plot(ridge, xvar = 'lambda')
#交叉验证
cv_ridge<-cv.glmnet(X,as.numeric(y), alpha=0)
cv_ridge$lambda.min
plot(cv_ridge)
#模型拟合
ridge_coef = coef(object = cv_ridge, s = cv_ridge$lambda.min)
ridge_coef



#############################################
#广义可加模型
#立方样条平滑cr+REML
############################################
#单因子广义可加模型
summary(gam(log(PM2.5)~s(PM10,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(SO2,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(NO2,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(CO,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(O3,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MAXT,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(AT,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MINT,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(AH,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MAXP,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(AP,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MINP,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MAXW,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(RAIN,bs="cr"),data=data_train,method="REML"))
summary(gam(log(PM2.5)~s(MAXR,bs="cr"),data=data_train,method="REML"))
# 剔除不显著的变量MAXW，观察到自由度均大于1，无线性部分

# 单因子图
plot(gam(log(PM2.5)~s(PM10,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(SO2,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(NO2,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(CO,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(O3,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MAXT,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(AT,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MINT,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(AH,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MAXP,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(AP,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MINP,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MAXW,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(RAIN,bs="cr"),data=data_train,method="REML"),shade=TRUE)
plot(gam(log(PM2.5)~s(MAXR,bs="cr"),data=data_train,method="REML"),shade=TRUE)

#多因素GAM模型
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(NO2,bs="cr")
           +s(CO,bs="cr")+s(O3,bs="cr")+s(MAXT,bs="cr")+s(AT,bs="cr")
           +s(MINT,bs="cr")+s(AH,bs="cr")+s(MAXP,bs="cr")+s(AP,bs="cr")
           +s(MINP,bs="cr")+s(RAIN,bs="cr")+s(MAXR,bs="cr"),data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
concurvity(gam_model,full=FALSE)#发现模型中很多变量不显著，存在共曲线性，要删减变量
#首先根据相关性，处理气温和气压部分，仅保留方差解释度最高的MINT
#利用筛选出的变量拟合多因素GAM模型
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(NO2,bs="cr")
           +s(CO,bs="cr")+s(O3,bs="cr")+s(MINT,bs="cr")+s(AH,bs="cr")           
           +s(RAIN,bs="cr")+s(MAXR,bs="cr"),data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
concurvity(gam_model,full=FALSE)
#共曲线性消失，但NO2 MINT AH、RAIN、MAXR自由度为1
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+NO2+s(CO,bs="cr")+s(O3,bs="cr")
           +MINT+AH+RAIN+MAXR,data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
#剔除不显著变量no2 mint rain maxr
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(CO,bs="cr")+s(O3,bs="cr")+AH,
              data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
#最终输入变量为：PM10 So2 Co O3 AH
plot(gam_model,shade=TRUE)#因素对PM2.5影响图

#预测
# 多元线性回归模型
pred1 <- predict(liner_step, newdata = data_test[,c('PM10','SO2','NO2','CO','O3','AP','RAIN','MAXR')])
pred1
ggplot(data = NULL, mapping = aes(pred1, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = '预测值', y = '实际值')
RMSE1 <- sqrt(mean((data_test$PM2.5-pred1)**2))
RMSE1
# BOX-COX多元线性回归模型
pred2 <- predict(liner_sqrt_step, newdata = data_test[,c('PM10','SO2','NO2','CO','O3','MINT','MAXW','RAIN','MAXR')])
pred2^2
ggplot(data = NULL, mapping = aes(pred2^2, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = '预测值', y = '实际值')
RMSE2 <- sqrt(mean((sqrt(data_test$PM2.5)-pred2)**2))
RMSE2
# 加权多元线性回归模型
pred3 <- predict(liner_weight_vif, newdata = data_test[,c('PM10','NO2','AH','MAXR')])
pred3
ggplot(data = NULL, mapping = aes(pred3, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = '预测值', y = '实际值')
RMSE3 <- sqrt(mean((data_test$PM2.5-pred3)**2))
RMSE3
# 岭回归模型
pred4 = predict(object = cv_ridge, s = cv_ridge$lambda.min, 
                newx = as.matrix(data_test[,3:17]))
ggplot(data = NULL, mapping = aes(pred4, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = '预测值', y = '实际值')
RMSE4 <- sqrt(mean((data_test$PM2.5-pred4)**2))
RMSE4
# 广义可加模型
pred5 = predict(gam_model,newdata = data_test[,3:17])
exp(pred5)
ggplot(data = NULL, mapping = aes(exp(pred5), data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = '预测值', y = '实际值')
RMSE5 <- sqrt(mean((log(data_test$PM2.5)-pred5)**2))
RMSE5






