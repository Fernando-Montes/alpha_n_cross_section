library(plotly)
# Parallel computing
library(doParallel)
registerDoParallel(cores=4)

# n1: average rate of room background event creating 1 neutron
# n2: average rate of room background event creating 2 neutrons
# b1: average rate of beam-induced background event creating 1 neutron
# b2: average rate of beam-induced background event creating 2 neutrons

do_experiment = function(n1, n2, b1, b2, e, l) {
  mult = c(0,0,0,0)
  for (i in 1:l) {
    n1_per_TW = rpois(1,n1) # Vector of l time-windows, number of n1 events created in each time-window
    n2_per_TW = rpois(1,n2) # Vector of l time-windows, number of n2 events created in each time-window
    b1_per_TW = rpois(1,b1) # Vector of l time-windows, number of b1 events created in each time-window
    b2_per_TW = rpois(1,b2) # Vector of l time-windows, number of b2 events created in each time-window
    # Vector of l time-windows, number of created neutrons in each time-window
    total_gen_neutrons = n1_per_TW + 2*n2_per_TW + b1_per_TW + 2*b2_per_TW 
    # Vector of l time-windows, number of detected neutrons in each time-window
    total_det_neutrons = rbinom( total_gen_neutrons, n = 1, prob = e )
    if ( total_det_neutrons == 0 ) mult[1] = mult[1] + 1
    else if ( total_det_neutrons == 1 ) mult[2] = mult[2] + 1
    else if ( total_det_neutrons == 2 ) mult[3] = mult[3] + 1
    else if ( total_det_neutrons == 3 ) mult[4] = mult[4] + 1
  }
  return( mult )
}

n1_mean = 0.008907
n2_mean = 623.303e-06

ptime <- system.time({
  table <- foreach(k=1:4, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
    return( do_experiment(n1_mean, n2_mean, 0, 0, 0.22, 89424000) )
  }
})[3]
print(ptime)

fileName <- "~/Dropbox/Courses/R/alpha_n/n1_n2.RData"
save(table, file = fileName)




tableTotal_exp = tableTotal[tableTotal$M1 < 505 & tableTotal$M1 > 500,]

p = ggplot(tableTotal_exp, aes(x=n1, y=n2)) + geom_point(size = 0.3) + 
  xlab("n1") + ylab("n2")
ggplotly( p + stat_density_2d() )

p = ggplot(data = data.frame(res = tableTotal$M1 ), aes(x = res)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30,
                 col = 'red', alpha = 0.5, fill = 'red') +
  geom_density(aes(y = ..density..), col = 'blue') + xlab('M1') 
ggplotly(p)
