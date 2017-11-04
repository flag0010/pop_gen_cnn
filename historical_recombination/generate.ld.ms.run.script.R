n.sims <- 210000

msSimulate  <- function(morgans.per.bp, ne, filename, n.chrom = 50, mu = 1.5e-8, bp = 20001){
  # converts parameters to hudson style, note nchrom, mu and bp are fixed
  total.genetic.distance 	<- morgans.per.bp * (bp - 1)
  rho 					<- 4 * ne * total.genetic.distance
  theta					<- 4 * ne * mu * bp
  k = sprintf("./ms %s 1 -t %s -r %s %s >> %s", n.chrom, theta, rho, bp, filename)
  return(k)
}


doSim 		<- function(index){
  # picks parameters
  mbp   <- 10^(runif(n = 1,min = -8, max = -6))
  my.ne <- sample(c(1000, 2000, 5000, 10000, 15000, 20000, 50000), 1)
  msSimulate(morgans.per.bp = mbp, ne = my.ne, filename = 'all.LD.sims.txt' )
}

x = unlist(lapply(1:n.sims, doSim))
write.table(x, 'run.ms.sh' sep='\n', quo=F, row=F)
