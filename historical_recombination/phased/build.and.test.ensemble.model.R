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
l = loess(score~nn,span=0.15)
#plot(real[order(real)], predict(l, real[order(real)]), ty='l', ylim=c(0,1))

#now extract weights for weighted mean from loess
#used the NN hear to guess what rho should be, because it had a greater accuracy 
#(it would be cheating to use the real rho value here, b/c in practice you'd never have that)
ensemble_weights = predict(l, nn)

#ensemble prediction using weighted mean
ensemble_pred = (nn*ensemble_weights + ldhat*(1-ensemble_weights))
#plot(real, ensemble_pred)

rmse_ensemble = mean((real-ensemble_pred)^2)^0.5
rmse_ldhat = mean((ldhat-real)^2)^0.5
rmse_nn = mean((nn-real)^2)^0.5

print (c('rmse ldhat', rmse_ldhat))
print (c('rmse ldhat', rmse_nn))
print (c('rmse ldhat', rmse_ensemble))

#[1] "rmse ldhat"         "0.0160305611954135"
#[1] "rmse ldhat"         "0.0125788914393887"
#[1] "rmse ldhat"         "0.0121620028651569"
