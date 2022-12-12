library(ggplot2)
library(reshape)
library(olsrr)
library(fastDummies)
library(rstan)
library(rstanarm)
library(SIBER)
library(dplyr)


set.seed(100)

# function built to load the data and do some transformations/create variables
load_data <- function(zipCodeDummy=FALSE){
  housingData <- read.csv("housingData.csv")
  housingData$date <- as.Date(housingData$date, "%Y%m%d")
  housingData$sales_year <- as.numeric(substring(housingData$date,1,4))
  housingData$age <- housingData$sales_year - housingData$yr_built
  housingData[housingData$yr_renovated==0, "was_renovated"] <- 0
  housingData[housingData$yr_renovated!=0, "was_renovated"] <- 1
  if(zipCodeDummy==TRUE){
    dummyzips <- dummy_cols(housingData$zipcode)
    housingData <- cbind(housingData, dummyzips[,c(2:ncol(dummyzips))])
    
  }
  housingData$logPrice <- log(housingData$price)
  
  housingData  
}

# function that creates a correlation matrix given data
create_correlation_matrix <- function(dataset){
  correlationMatrix <- data.frame()
  for(indivCol in c(3:ncol(dataset))){
    for(indivCol2 in c(3:ncol(dataset))){
      factor1 <- names(dataset)[indivCol]
      factor2 <- names(dataset)[indivCol2]
      correlation <- cor(dataset[,indivCol], dataset[,indivCol2], use="everything")
      temp <- data.frame(factor1, factor2, correlation)
      print(temp)
      correlationMatrix <- rbind(correlationMatrix, temp)
      
    }
  }
  correlationMatrix
}
# Runs variable selection algorithm
run_linear_model_variable_Selection <- function(dataset){
  linearModel <- lm(logPrice~age + bathrooms + bedrooms + condition + floors +
                      grade + sqft_above + sqft_basement + sqft_living + sqft_lot +
                      view + was_renovated + waterfront, data = housingData)
  k <- ols_step_all_possible(linearModel)
  k
}

# creates data exploration charts
plot_data_exploration_charts <- function(dataset){
  meltedData <- melt(dataset, id=c("id", "date", "price"))
  pdf("dataExploration.pdf", height=8.5, width=11)
  a <- ggplot(meltedData, aes(x=value, y=price))+
    facet_wrap(~variable, nrow=5, ncol=5, scales="free_x")+geom_point()+geom_smooth(method='lm')
  print(a)
  dev.off()

}
# normalizes the non dummy-variable values
normalizeHousingData <- function(dataset){
  numericalFields <- c("bedrooms", "bathrooms", "sqft_living", "sqft_lot", 
                       "floors", "view", "condition", "grade", "sqft_above",
                       "sqft_basement", "sqft_living15", "sqft_lot15", "age")
  
  dataset[, numericalFields] <- scale(dataset[, numericalFields])
  dataset
}



housingData <- load_data(zipCodeDummy=TRUE)
housingData <- normalizeHousingData(housingData)

# 1-Factor Linear Fit
linearModel1Factor <- lm(logPrice~sqft_living, data=housingData)
# Scatterplot of the fit
plot(housingData$logPrice~housingData$sqft_living, main="Scatter Plot and Linear Fit",
     xlab="sqft_living", ylab="log(Home Price")
abline(a=linearModel1Factor$coefficients[1],
       b=linearModel1Factor$coefficients[2], col="red")






# creates PDF of diagnostic fit plots
pdf("linearModel1FactorCharts.pdf")
plot(linearModel1Factor)
dev.off()



# Bayesian 1-Factor Regression
bayesModel1Factor <- stan_glm(logPrice ~ sqft_living, data=housingData)

# Calculates the Bayesian R-Sqaured for the 1-factor fit.
bayesR2.1Factor <- bayes_R2(bayesModel1Factor)

#bayes R2 decomposition
y_pred <- posterior_linpred(bayesModel1Factor)
var_fit <- apply(y_pred,1,var)
var_res <- as.matrix(bayesModel1Factor, pars=c("sigma"))^2

hist(var_fit, main="Explained Variance Statistic histogram",
     xlab="Variance(E((pred(y_n)|θ)))")

hist(var_res, main="Residual Variance Statistic histogram",
     xlab="Variance(E((y_n-pred(y_n)|θ)))")


# plot histogram of R-2 from Posterior Draws
hist(bayesR2.1Factor,
     main="Bayesian R2 for 1-Factor Model Price~sqft_living",
     xlab="R-Squared"
)

# Trace plots to show convergence of the algorithm
plot(bayesModel1Factor, "trace")

bayesModel1Factor.posterior <- as.matrix(bayesModel1Factor)

mcmc_hist(bayesModel1Factor.posterior,pars=c("sqft_living", "sigma"), prob=0.95) + ggtitle("Posterior Distributions",
                                                             "with medians and 95% intervals")

# code to draw scatter plot and posterior regression lines
n_draws<-2000
alpha_level <- .15
color_draw <- "grey60"
color_mean <- "red"
bayesModel1FactorDraws <- bayesModel1Factor %>% as_tibble()
pdf("bayesianScatterPlot.pdf", height=8.5, width=11)
a <- ggplot(housingData, aes(x=sqft_living, y=logPrice)) + 
  geom_abline(aes(intercept=`(Intercept)`, slope=sqft_living),
              data=sample_n(bayesModel1FactorDraws, n_draws),
              color=color_draw,
              alpha=alpha_level) + 
  geom_abline(intercept=mean(bayesModel1FactorDraws$`(Intercept)`), 
              slope=mean(bayesModel1FactorDraws$sqft_living),
              size=1,
              color=color_mean) + geom_point() + 
  labs(x = "standardized sqft living space",
       y = "log(home prices)",
       title="Visualization of 2000 Regression Lines from the Posterior Distribution")
print(a)
dev.off()

# Posterior distribution mean vs data mean
pp_check(bayesModel1Factor, "stat") + 
  ggtitle("Comparing Mean Log(Home Prices) from Posterior Draws versus the Data") + 
  labs(x="Log(Home Prices)")

# Posterior density overlay plot
pp_check(bayesModel1Factor, "dens_overlay") + 
  ggtitle("Comparing Posterior Distribution to Data Density")

#Posterior draws 2-D plot
pp_check(bayesModel1Factor, "stat_2d") + ggtitle("Comparing mean and Standard 
                                                 Deviaton of replicated versus 
                                                 observed data")


linearModelBest <- lm(logPrice~age + sqft_living + grade ++ .data_98001+
                        .data_98002+.data_98003+.data_98004+.data_98005+
                        .data_98006+.data_98007+.data_98008+.data_98010+
                        .data_98011+.data_98014+.data_98019+.data_98022+
                        .data_98023+.data_98024+.data_98027+.data_98028+
                        .data_98029+.data_98030+.data_98031+.data_98032+
                        .data_98033+.data_98034+.data_98038+.data_98039+
                        .data_98040+.data_98042+.data_98045+.data_98052+
                        .data_98053+.data_98055+.data_98056+.data_98058+
                        .data_98059+.data_98065+.data_98070+.data_98072+
                        .data_98074+.data_98075+.data_98077+.data_98092+
                        .data_98102+.data_98103+.data_98105+.data_98106+
                        .data_98107+.data_98108+.data_98109+.data_98112+
                        .data_98115+.data_98116+.data_98117+.data_98118+
                        .data_98119+.data_98122+.data_98125+.data_98126+
                        .data_98133+.data_98136+.data_98144+.data_98146+
                        .data_98148+.data_98155+.data_98166+.data_98168+
                        .data_98177+.data_98178+.data_98188+.data_98198+
                        .data_98199, data= housingData)
bayesModelBest <- stan_glm(logPrice~age + sqft_living + grade ++ .data_98001+
                             .data_98002+.data_98003+.data_98004+.data_98005+
                             .data_98006+.data_98007+.data_98008+.data_98010+
                             .data_98011+.data_98014+.data_98019+.data_98022+
                             .data_98023+.data_98024+.data_98027+.data_98028+
                             .data_98029+.data_98030+.data_98031+.data_98032+
                             .data_98033+.data_98034+.data_98038+.data_98039+
                             .data_98040+.data_98042+.data_98045+.data_98052+
                             .data_98053+.data_98055+.data_98056+.data_98058+
                             .data_98059+.data_98065+.data_98070+.data_98072+
                             .data_98074+.data_98075+.data_98077+.data_98092+
                             .data_98102+.data_98103+.data_98105+.data_98106+
                             .data_98107+.data_98108+.data_98109+.data_98112+
                             .data_98115+.data_98116+.data_98117+.data_98118+
                             .data_98119+.data_98122+.data_98125+.data_98126+
                             .data_98133+.data_98136+.data_98144+.data_98146+
                             .data_98148+.data_98155+.data_98166+.data_98168+
                             .data_98177+.data_98178+.data_98188+.data_98198+
                             .data_98199, data= housingData, iter=100000, chains=6)
# Multi-factor Bayesian R-Squared 
y_pred_best <- posterior_linpred(bayesModelBest, draws=4000)
var_fit_best <- apply(y_pred,1,var)
var_res_best <- (as.matrix(bayesModelBest, pars=c("sigma"))[1:4000])^2
bayesR2.best <- var_fit_best / (var_fit_best+var_res_best)

hist(bayesR2.best,
     main="Bayesian R2 for 4-Factor Model Price~sqft_living + age + grade + zipcode dummy",
     xlab="R-Squared"
)

# Posterior distribution mean vs data mean
pp_check(bayesModelBest, "stat", ndraws=100) + 
  ggtitle("Comparing Mean Log(Home Prices) from Posterior Draws versus the Data") + 
  labs(x="Log(Home Prices)")


# Posterior density overlay plot
pp_check(bayesModelBest, "dens_overlay", ndraws=5) + 
  +     ggtitle("Comparing Mean Log(Home Prices) from Posterior Draws versus the Data") + 
  +     labs(x="Log(Home Prices)")


pdf("linearModelBestCharts.pdf")
plot(linearModelBest)
dev.off()
