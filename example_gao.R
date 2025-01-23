small_data <- readRDS("C:/Users/WD/Downloads/MiGA_eQTL.chr2_ENSG00000151694.univariate_data.rds")


y= small_data$ENSG00000151694$residual_Y[[3]]
X= small_data$ENSG00000151694$residual_X[[3]]
res_susie =susieR::susie(X = X,y=y,L=20   )
res_susie$sets

plot(y , predict(res_susie , X) )

plot( res_susie$pip)
res_susie$sigma2/var(y)
library(susieRsmall)

  res_susie_small =susieRsmall::susie(X = X,y=y,L=10, max_iter = 100 , estimate_prior_method = "EM")
res_susie_small$sets

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
res_susie_small$mu[1:10,1:10]
res_susie  =susieR ::susie(X = X[,1:1000],y=y,L=10, max_iter =20)
res_susie $sets

res_susie $mu[1:10,1:10]


lm( y~X[,1])

plot(res_susie $pip)

points(res_susie_small$pip, pch=20, col=3)



res_susie$sigma2
res_susie_small$sigma2

var(X[,1])
var(y)
res_susie$V
res_susie_small$V

o=o+1
