

#1. Read input
input <- read.csv(".\\input.csv", sep = ";")
input <- input[,c(5,9,10,19,20,22,24:27)]
input$Kbin_1 <- as.factor(input$Kbin_1)

#2. Model fit
#partition dataset into training and testing (70 vs 30%)
library(caret)
Train <- createDataPartition(input$rep, p=0.7, list=FALSE)
training <- input[ Train, ]
testing <- input[ -Train, ]
#fit model
mod_fit <- glm(rep ~ bio01_1 + forestCo_1 + slope_1 + soil_1 + vpd_1 + Kbin_1,  data=training, family="binomial")

#3. Basic model eval
#check variable importance
summary(mod_fit) #is variable contribution significant?
#Binned residual plots
library(arm)
binnedplot(fitted(mod_fit), 
           residuals(mod_fit, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")

#plotting response curves (repeat for each significant variable)
glm.probs <- predict(mod_fit,newdata=testing, type = "response") #Use model to predict rep mode on testing dataset
hist(glm.probs[testing$rep == "Larv"])
hist(glm.probs[testing$rep == "Viv"])

training$prob <- fitted(mod_fit)
training$type_num <- ifelse(training$rep == "Viv", 1, 0)

library(ggplot2)
ggplot(training, aes(bio01_1, type_num)) +
  geom_point() +
  stat_smooth(aes(bio01_1, prob), se = F) +
  labs(x = "bio01", y = "Pr(Pueriparous)", title = "larviparous vs. Pueriparous ~ Average Mean Temperature (bio01)")
#for inclusion of Kbin in response curve plots
ggplot(training, aes(soil_1, type_num)) +
  geom_point() +
  stat_smooth(aes(soil_1, prob, col = Kbin_1), se = F) +
  labs(x = "soil", y = "Pr(Pueriparous)", title = "Larviparous vs. Pueriparous ~ Soil moisture vs. Lithology")

##4. In depth model evaluation
#4.1 Discrimination
#4.1.1 tss (not kept for publication)
#thresh <- seq(0.01,1,0.01)
#id <- c("Viv", "Larv")
#tss.results <- data.frame("threshold" = rep(0, length(thresh)), "sensitivity" = rep(0, length(thresh)), "specificity" = rep(0, length(thresh)), "maxTSS" = rep(0, length(thresh)))
#for (i in 1: length(thresh))
#     {
#  glm.pred <- ifelse(glm.probs > thresh[i], id[1], id[2])
#  cm <- table(glm.pred, testing$mode) #Confusion matrix
#  if (dim(cm)[1] == 1)
#  {
#    cm2 <- c(0,0)
#    cm <- rbind(cm, cm2)
#    rm(cm2)
#  }
#  sp <- cm[1,1]/(cm[1,1] + cm[2,1])
#  sn <- cm[2,2]/(cm[1,2] + cm[2,2])
#  tss <- sn + sp -1
#  tss.results[i,] <- c(thresh[i], sn, sp, tss)
#  rm(sp, sn, tss)
#}
#tss.results[which(tss.results$maxTSS == max(tss.results$maxTSS)),]
#Create confusion matrix for max TSS
#glm.pred.best <- ifelse(glm.probs > thresh[which(tss.results$maxTSS == max(tss.results$maxTSS))], id[1], id[2])
#cm.best <- table(glm.pred.best, testing$rep) #Confusion matrix
#cm.best
#rm(cm, sn, sp, tss, tss.results, glm.pred, glm.pred.best,cm.best)

#4.1.2 auc
#for all
library(ROCR)
#test data
prob <- predict(mod_fit, newdata=testing, type="response")
pred <- prediction(prob, testing$rep)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc.test <- auc@y.values[[1]]
auc.test
#training data
prob <- predict(mod_fit, newdata=training, type="response")
pred <- prediction(prob, training$rep)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc.training <- auc@y.values[[1]]
auc.training

#4.2 Calibration
library(pscl)
pR2(mod_fit)  # Pseudo R^2.look for 'McFadden'. Closer to 1 is better

#5 Save final model
training.Larviv <- training
testing.Larviv <- testing
mod_fit.Larviv <- mod_fit


###subgroups (repeat for each pairwise model)
#1. Read input
input <- read.csv(".\\input.csv", sep = ";")
input <- input[,c(5,9,10,19,20,22,24:27)]
input$Kbin_1 <- as.factor(input$Kbin_1)

#2. Model fit
input.sub <- subset(input, input$mode == "Viv_Salg" | input$mode == "Viv_Lycia")
input.sub$mode <- factor(input.sub$mode)
#partition dataset into training and testing (70 vs 30%)
library(caret)
Train <- createDataPartition(input.sub$mode, p=0.7, list=FALSE)
training <- input.sub[ Train, ]
testing <- input.sub[ -Train, ]
#fit model (for stepwise models, manually edit to include only variables of interest)
mod_fit <- glm(mode ~  Kbin_1 + forestCo_1 + slope_1 + bio01_1 + soil_1,  data=training, family="binomial", maxit = 100)

#3. Basic model eval
#check variable importance
summary(mod_fit) #is value contribution significant?
#check residuals. 95% of points should fall within confidence interval (the grey lines)
library(arm)
binnedplot(fitted(mod_fit), 
           residuals(mod_fit, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")
#plotting response curves (repeat for each significant variable)
glm.probs <- predict(mod_fit,newdata=testing, type = "response") #Use model to predict rep mode on testing dataset
hist(glm.probs[testing$mode == "Viv_Lycia"])
hist(glm.probs[testing$mode == "Viv_Salg"])

training$prob <- fitted(mod_fit)
training$type_num <- ifelse(training$mode == "Viv_Lycia", 1, 0)

library(ggplot2)
ggplot(training, aes(soil_1, type_num)) +
  geom_point() +
  stat_smooth(aes(soil_1, prob), se = F) +
  labs(x = "soil", y = "Pr(Lyciasalamandra)", title = "Alpine salamanders vs. Lyciasalamandra ~ Soil moisture (soil)")

#for inclusion of Kbin in response curve plots
ggplot(training, aes(soil_1, type_num)) +
  geom_point() +
  stat_smooth(aes(soil_1, prob, col = Kbin_1), se = F) +
  labs(x = "soil", y = "Pr(S. algira)", title = "Lyciasalamandra vs. Salamandra algira ~ Soil moisture (soil) vs. lithology")

##4. In depth model evaluation
#4.1 Discrimination
#4.1.1 tss (not kept for publication)
#thresh <- seq(0.01,0.99,0.01)
#id <- c("Viv_Lycia", "Viv_Salg")
#tss.results <- data.frame("threshold" = rep(0, length(thresh)), "sensitivity" = rep(0, length(thresh)), "specificity" = rep(0, length(thresh)), "maxTSS" = rep(0, length(thresh)))
#for (i in 1: length(thresh))
#{
#  glm.pred <- ifelse(glm.probs > thresh[i], id[1], id[2])
#  cm <- table(glm.pred, testing$mode) #Confusion matrix
#  if (dim(cm)[1] == 1)
#  {
#    cm2 <- c(0,0)
#    cm <- rbind(cm, cm2)
#    rm(cm2)
#  }
#  sp <- cm[1,1]/(cm[1,1] + cm[2,1])
#  sn <- cm[2,2]/(cm[1,2] + cm[2,2])
#  tss <- sn + sp -1
#  tss.results[i,] <- c(thresh[i], sn, sp, tss)
#  rm(sp, sn, tss)

#}
#tss.results[which(tss.results$maxTSS == max(tss.results$maxTSS)),]
#glm.pred.best <- ifelse(glm.probs > thresh[which(tss.results$maxTSS == max(tss.results$maxTSS))[1]], id[1], id[2])
#cm.best <- table(glm.pred.best, testing$mode) #Confusion matrix
#cm.best
#rm(cm, sn, sp, tss, tss.results, glm.pred, glm.pred.best,cm.best)
#4.1.2 auc
#for all
library(ROCR)
#test data
prob <- predict(mod_fit, newdata=testing, type="response")
pred <- prediction(prob, testing$mode)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc.test <- auc@y.values[[1]]
auc.test

#training data
prob <- predict(mod_fit, newdata=training, type="response")
pred <- prediction(prob, training$mode)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
auc <- performance(pred, measure = "auc")
auc.training <- auc@y.values[[1]]
auc.training

#4.2 Calibration
library(pscl)
pR2(mod_fit)  # Pseudo R^2.look for 'McFadden'. Closer to 1 is better


#5 Save final model #edit object names mannually for each model iteration you wish to save
training.LyciaSalg <- training
testing.LyciaSalg <- testing
mod_fit.LyciaSalg <- mod_fit
