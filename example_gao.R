rm(list=ls())

small_data <- readRDS("C:/Users/WD/Downloads/MiGA_eQTL.chr2_ENSG00000151694.univariate_data.rds")


y= small_data$ENSG00000151694$residual_Y[[3]]
X= small_data$ENSG00000151694$residual_X[[3]]
res_susie =susieR::susie(X = X,y=y,L=20   )
res_susie$sets

plot(y , predict(res_susie , X) )

plot( res_susie$pip)
res_susie$sigma2/var(y)
library(susieRsmall)

res_susie_small =susieRsmall::susie(X = X,y=y,L=10, max_iter = 20 , estimate_prior_method = "EM")

res_susie_small$V
res_susie_small$sets
susie_get_cs(res_susie_small,
             X = X,
             coverage = 0.9,
             min_abs_corr = 0.1)



plot( res_susie$pip,
      res_susie_small$pip)


plot(predict(res_susie_small,X),y)
points(predict(res_susie,X),y)
plot(res_susie_small$mu[1,],res_susie$mu[1,])
abline(a=0,b=1
       )


plot( exp(res_susie_small$lbf_variable[1,])/sum(exp(res_susie_small$lbf_variable[1,])))

cor(X[ ,which(res_susie_small$alpha[1,]>0.05)])^2

plot(apply( res_susie_small$alpha,2, mean))

h2_Gaussian = 100*(1-res_susie$sigma2/var(y))

h2_SS = 100*(1-res_susie_small$sigma2/var(y))
par(mfrow=c(1,2))
susie_plot( res_susie, y="PIP",
            main=paste( "SuSiE Gaussian SER\n
            CSs explain",
                        round( c(h2_Gaussian),digits=2), "% of the variance"))

susie_plot( res_susie_small, y="PIP",
            main=paste( "SuSiE Servin Stephens SER\n
            CS explain",
                        round( c(h2_SS),digits=2), "% of the variance"))



par(mfrow=c(1,1))

res_susie$sets[[1]][[1]]
res_susie_small$lbf
res_susie_small$sets[[1]][[1]]
plot(res_susie_small$alpha[1,])

res_susie$V
res_susie_small$V

plot(res_susie $pip)

points(res_susie_small$pip, pch=20, col=3)
res_susie_small$sets
res_susie$elbo
res_susie_small$elbo

var(predict(res_susie_small,X),y)

res_susie_small$sigma2

res_susie_small$sigma2

var(y)
res_susie_small$lbf
res_susie $lbf

plot(y , predict(res_susie_small, X) )


res_susie_small$V

res_list= list()

for ( i in 1:ncol(X)){


  res_list[[i]]= summary(lm(y~X[,i]))$coefficients [2,]
  print( i)
}
marginal = do.call(rbind, res_list)
plot( -log10 (marginal[,4]))





#### too look at ---



o=10

set.seed(o)
y = X[,1]+ rnorm( nrow(X),sd= sd(X[,1]))


devtools::load_all()
res_susie_small =susieRsmall::susie(X = X[,1:1000],y=y,L=10,
                                    max_iter =20)
res_susie_small =susieRsmall::susie(X = X[,1:1000],y=y,L=10,
                                    max_iter =20,
                                    estimate_prior_method = "EM")
res_susie_small$sets
res_susie  =susieR ::susie(X = X[,1:1000],y=y,L=10, max_iter =20)
res_susie $sets



plot(res_susie $pip)

points(res_susie_small$pip, pch=20, col=3)



res_susie$sigma2
res_susie_small$sigma2

var(X[,1])
var(y)
res_susie$V
res_susie_small$V

o=o+1

