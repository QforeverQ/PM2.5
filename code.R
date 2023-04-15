library(ggplot2)
library(gridExtra)
library(glmnet)
library(car)   # ģ�����
library(mgcv)  # ����ɼ�ģ��

setwd("D:/�о���/ͳ�ƻ���/��ĩ����ҵ")
# ����ѵ�����Ͳ��Լ�����
data_train <- read.csv("D:/�о���/ͳ�ƻ���/��ĩ����ҵ/����pm2.5����ѵ����.csv",header=TRUE,fileEncoding='UTF-8')
data_test <- read.csv("D:/�о���/ͳ�ƻ���/��ĩ����ҵ/����pm2.5���ݲ��Լ�.csv",header=TRUE)
#�����ݱ�Ϊ���ݿ�
data_train<-as.data.frame(data_train)
data_test<-as.data.frame(data_test)

head(data_train)
attach(data_train)
# �鿴���ݽṹ
str(data_train)
# �鿴���ݵĸ�����Ϣ
summary(data_train)

# PM2.5��15��Э�������ɢ��ͼ
for (i in 1:15) {
  pg <- ggplot(data = data_train, mapping = aes(x = data_train[,i+2], y = data_train[,2])) + 
    geom_point(fill = 'steelblue') + labs(x = names(data_train)[i+2], y = 'PM2.5')
  ggsave(paste0(i, "_", ".png"), pg, width = 3, height = 3)
}

# ��ģ����Ԫ���Իع�ģ��
liner = lm(PM2.5 ~ PM10 + SO2 + NO2 + CO + O3 + MAXT + AT + MINT + AH + MAXP + AP +
             MINP + MAXW + RAIN + MAXR, data = data_train)
summary(liner)

# �𲽻ع���ɱ���ѡ��
liner_step <- step(liner)
summary(liner_step)


# ģ�����
# ��̬�Լ���
# ͼ�η�
#����ֱ��ͼ
hist(x = data_train$PM2.5, freq = FALSE, main = 'PM2.5��ֱ��ͼ',
     ylab = '���ܶ�ֵ',xlab = NULL, ylim = c(0,0.02), col = 'steelblue')
#���Ӻ��ܶ�ͼ
lines(density(data_train$PM2.5), col = 'red', lty = 1, lwd = 2)
#������̬�ֲ�ͼ
x <- data_train$PM2.5[order(data_train$PM2.5)]
lines(x, dnorm(x, mean(x), sd(x)),
      col = 'black', lty = 2, lwd = 2.5)

#����ͼ��
legend('topright',legend = c('���ܶ�����','��̬�ֲ�����'),
       col = c('red','black'), lty = c(1,2),
       lwd = c(2,2.5), bty = 'n')

# PPͼ
real_dist <- ppoints(data_train$PM2.5)
theory_dist <- pnorm(data_train$PM2.5, mean = mean(data_train$PM2.5), 
                     sd = sd(data_train$PM2.5))

# ��ͼ
plot(sort(theory_dist), real_dist, col = 'steelblue', 
     pch = 20, main = 'PPͼ', xlab = '�����ۼƸ���', 
     ylab = 'ʵ���ۼƸ���')

# ���ӶԽ�����Ϊ�ο���
abline(a = 0,b = 1, col = 'red', lwd = 2)

# QQͼ
qqnorm(data_train$PM2.5, col = 'steelblue', pch = 20,
       main = 'QQͼ', xlab = '���۷�λ��', 
       ylab = 'ʵ�ʷ�λ��')
# ���Ʋο���
qqline(data_train$PM2.5, col = 'red', lwd = 2)

# ͳ�Ʒ�
shapiro.test(data_train$PM2.5)

# ��PM2.5��������̬�ֲ�ʱ
# box-cox�任
powerTransform(liner_step)   
#������
liner_sqrt <- lm(sqrt(PM2.5) ~ PM10 + SO2 + NO2 + CO + O3 + MAXT + AT + MINT + AH + MAXP + AP +
                  MINP + MAXW + RAIN + MAXR, data = data_train)
summary(liner_sqrt)
liner_sqrt_step <- step(liner_sqrt)
summary(liner_sqrt_step)

# ���ع����Լ���
vif(liner_sqrt_step)  # ͨ��
plot(liner_sqrt_step)



# �������Լ��飺��ͨ��
ncvTest(liner_sqrt_step)

# ��Ȩ��С���˷�
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
powerTransform(liner_weight_step)  #����ҪBOX-COX�任
vif(liner_weight_step)  #�й�������

# ���ݱ��������Ժ͹������Խ��б���ɸѡ
liner_weight_vif <- lm(PM2.5 ~ PM10 + NO2 + AH + MAXR, data = data_train_W)
summary(liner_weight_vif)
vif(liner_weight_vif)

plot(liner_weight_vif)



#����ת��Ϊ����
M = as.matrix(data_train)
X = M[,3:17]
X=apply(X,2,as.numeric)
y = M[,2]

#��ع飺Ч����
#���ӻ�ȷ��lambda
ridge<-glmnet(X,y,alpha=0)
coef(ridge)
plot(ridge, xvar = 'lambda')
#������֤
cv_ridge<-cv.glmnet(X,as.numeric(y), alpha=0)
cv_ridge$lambda.min
plot(cv_ridge)
#ģ�����
ridge_coef = coef(object = cv_ridge, s = cv_ridge$lambda.min)
ridge_coef



#############################################
#����ɼ�ģ��
#��������ƽ��cr+REML
############################################
#�����ӹ���ɼ�ģ��
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
# �޳��������ı���MAXW���۲쵽���ɶȾ�����1�������Բ���

# ������ͼ
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

#������GAMģ��
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(NO2,bs="cr")
           +s(CO,bs="cr")+s(O3,bs="cr")+s(MAXT,bs="cr")+s(AT,bs="cr")
           +s(MINT,bs="cr")+s(AH,bs="cr")+s(MAXP,bs="cr")+s(AP,bs="cr")
           +s(MINP,bs="cr")+s(RAIN,bs="cr")+s(MAXR,bs="cr"),data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
concurvity(gam_model,full=FALSE)#����ģ���кܶ���������������ڹ������ԣ�Ҫɾ������
#���ȸ�������ԣ��������º���ѹ���֣�������������Ͷ���ߵ�MINT
#����ɸѡ���ı�����϶�����GAMģ��
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(NO2,bs="cr")
           +s(CO,bs="cr")+s(O3,bs="cr")+s(MINT,bs="cr")+s(AH,bs="cr")           
           +s(RAIN,bs="cr")+s(MAXR,bs="cr"),data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
concurvity(gam_model,full=FALSE)
#����������ʧ����NO2 MINT AH��RAIN��MAXR���ɶ�Ϊ1
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+NO2+s(CO,bs="cr")+s(O3,bs="cr")
           +MINT+AH+RAIN+MAXR,data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
#�޳�����������no2 mint rain maxr
gam_model=gam(log(PM2.5)~s(PM10,bs="cr")+s(SO2,bs="cr")+s(CO,bs="cr")+s(O3,bs="cr")+AH,
              data=data_train,method="REML")
summary(gam_model)
concurvity(gam_model)
#�����������Ϊ��PM10 So2 Co O3 AH
plot(gam_model,shade=TRUE)#���ض�PM2.5Ӱ��ͼ

#Ԥ��
# ��Ԫ���Իع�ģ��
pred1 <- predict(liner_step, newdata = data_test[,c('PM10','SO2','NO2','CO','O3','AP','RAIN','MAXR')])
pred1
ggplot(data = NULL, mapping = aes(pred1, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = 'Ԥ��ֵ', y = 'ʵ��ֵ')
RMSE1 <- sqrt(mean((data_test$PM2.5-pred1)**2))
RMSE1
# BOX-COX��Ԫ���Իع�ģ��
pred2 <- predict(liner_sqrt_step, newdata = data_test[,c('PM10','SO2','NO2','CO','O3','MINT','MAXW','RAIN','MAXR')])
pred2^2
ggplot(data = NULL, mapping = aes(pred2^2, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = 'Ԥ��ֵ', y = 'ʵ��ֵ')
RMSE2 <- sqrt(mean((sqrt(data_test$PM2.5)-pred2)**2))
RMSE2
# ��Ȩ��Ԫ���Իع�ģ��
pred3 <- predict(liner_weight_vif, newdata = data_test[,c('PM10','NO2','AH','MAXR')])
pred3
ggplot(data = NULL, mapping = aes(pred3, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = 'Ԥ��ֵ', y = 'ʵ��ֵ')
RMSE3 <- sqrt(mean((data_test$PM2.5-pred3)**2))
RMSE3
# ��ع�ģ��
pred4 = predict(object = cv_ridge, s = cv_ridge$lambda.min, 
                newx = as.matrix(data_test[,3:17]))
ggplot(data = NULL, mapping = aes(pred4, data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = 'Ԥ��ֵ', y = 'ʵ��ֵ')
RMSE4 <- sqrt(mean((data_test$PM2.5-pred4)**2))
RMSE4
# ����ɼ�ģ��
pred5 = predict(gam_model,newdata = data_test[,3:17])
exp(pred5)
ggplot(data = NULL, mapping = aes(exp(pred5), data_test$PM2.5)) + 
  geom_point(color = 'red', shape = 19) + 
  geom_abline(slope = 1, intercept = 0, size = 1) +
  labs(x = 'Ԥ��ֵ', y = 'ʵ��ֵ')
RMSE5 <- sqrt(mean((log(data_test$PM2.5)-pred5)**2))
RMSE5





