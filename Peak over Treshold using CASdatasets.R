library(CASdatasets)
library(MASS)

# Aim of this program is to evaluate the expected value of claims greater than chosen treshold so E(Y|Y>u)

##Loading data set

data(swmotorcycle)
data_set=swmotorcycle
str(data_set)
data_set<-data_set[which(data_set$Exposure>0),]
data_set=data_set[data_set$ClaimAmount>0,]
dim(data_set)


library(POT)
library(tea)
mrlplot(data_set$ClaimAmount)
abline(40000,0)
tcplot(data_set$ClaimAmount)

gg_plot_object <- ggplot(data_set$ClaimAmount, nexceed = min(data_set$ClaimAmount) - 1)
threshold=gg_plot_object$threshold #chosen threshold for distribution estimation 

fitted_gdp = fitgpd(data_set$ClaimAmount, threshold=threshold)
plot(fitted_gdp)
qq(fitted_gdp)
pp(fitted_gdp)
dens(fitted_gdp)
confint(fitted_gdp, prob = 0.95)

tau = fitted_gdp$param[[1]] 
ksi = fitted_gdp$param[[2]]

n = length(data_set$ClaimAmount)
n_threshold = length(data_set[data_set$ClaimAmount>threshold,"ClaimAmount"]) 

# we set chosen threshold to 80k
u = 80000
p=1-(1-pgpd(threshold, scale=tau, shape=ksi))*(1+ksi*(u-threshold)/tau)^(-1/ksi)

pgpd(u, loc=threshold, scale=tau, shape=ksi)
pgpd(threshold, loc=threshold, scale=tau, shape=ksi)

quantile_fitted_p = threshold+tau/ksi*((n*(1-p)/n_threshold)^(-1*ksi)-1)

(quantile_fitted_p+tau-ksi*threshold)/(1-ksi)

plot(1:100000,dgpd(1:100000,scale=tau, shape=ksi))
summary(data_set$ClaimAmount)

Y = data_set$ClaimAmount
format(Y, scientific = FALSE)
bins = seq(min(Y), max(Y), length.out = 11)
#hist(Y, breaks = bins, labels=bins, probability =TRUE, labels=bins)
plot(ecdf(Y[Y>threshold]))
points(Y[Y>threshold], pgpd(Y[Y>threshold], loc=threshold, scale=tau, shape=ksi), col="blue")
abline(v=threshold)
points(Y[Y>u], pgpd(Y[Y>u], loc=threshold, scale=tau, shape=ksi), col="green")
abline(v=u, col="green") #chosen threshold for data


mean(rgpd(1000000, loc=u, scale=tau, shape=ksi)) # estimated E(Y|Y>u)
e_u = (u+tau-ksi)/(1-ksi) # Real E(Y|Y>u)
(mean(Y[Y>u])-e_u)/mean(Y[Y>u])
tau/(1-ksi)
abline(v=e_u, col="red") #chosen threshold for data

