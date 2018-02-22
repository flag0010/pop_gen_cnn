#read in data
x = read.csv('input.data.for.ensemble.csv')
real= x[,1]
nn = x[,2]
ldhat = x[,3]

#build scoring func 
#give 1 to cases where nn better than ldhat, and zero for reverse
ldhat_abs_resid = abs(ldhat-real)
nn_abs_resid = abs(nn-real)
score = as.numeric(ldhat_abs_resid > nn_abs_resid)

#fit a loess func (this is a little different than the version in python, but same idea)
#used the NN here and below to guess what rho should be, because it had a greater accuracy 
#(it would be cheating to use the real rho value here, b/c in practice you'd never have that outside of simulations)
l = loess(score~nn,span=0.15)
#plot(real[order(real)], predict(l, real[order(real)]), ty='l', ylim=c(0,1))

#now extract weights for weighted mean from loess
ensemble_weights = predict(l, nn)

#ensemble prediction using weighted mean
ensemble_pred = (nn*ensemble_weights + ldhat*(1-ensemble_weights))
#plot(real, ensemble_pred)

rmse_ensemble = mean((real-ensemble_pred)^2)^0.5
rmse_ldhat = mean((ldhat-real)^2)^0.5
rmse_nn = mean((nn-real)^2)^0.5

r2_ensemble = cor(real, ensemble_pred)^2
r2_ldhat = cor(real, ldhat)^2
r2_nn = cor(real, nn)^2

print (c('rmse ldhat', rmse_ldhat))
print (c('rmse NN', rmse_nn))
print (c('rmse ENSEMBLE', rmse_ensemble))

print (c('r2 ldhat', r2_ldhat))
print (c('r2 NN', r2_nn))
print (c('r2 ENSEMBLE', r2_ensemble))

#[1] "rmse ldhat"         "0.0160305611954135"
#[1] "rmse NN"         "0.0125788914393887"
#[1] "rmse ENSEMBLE"         "0.0121620028651569"
#[1] "r2 ldhat"          "0.768138741795386"
#[1] "r2 NN"             "0.861572346243779"
#[1] "r2 ENSEMBLE"       "0.881294110333308"
