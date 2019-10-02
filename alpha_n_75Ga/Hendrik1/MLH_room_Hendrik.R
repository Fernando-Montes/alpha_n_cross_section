# Combinatorial function
comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}

# Function returning probability of detecting det number of neutrons in time-window coming from s1
prob_s1 = function(det, s1, eff) {
  prob = 0
  for(i in det:5) {
    prob = prob + dpois(i, s1)*eff^det*(1-eff)^(i-det)*comb(i, det)
  }
  return( prob )
}

# Function returning probability of detecting det number of neutrons in time-window coming from s2
prob_s2 = function(det, s2, eff) {
  prob = 0
  for(i in ceiling(det/2):5) {
    prob = prob + dpois(i, s2)*eff^det*(1-eff)^(2*i-det)*comb(2*i, det)
  }
  return( prob )
}

# Function returning probability of 0 neutron detected in time-window
prob_det0 = function(n1, n2, eff){
  return( prob_s1(0, n1, eff)*prob_s2(0, n2, eff) )
}

# Function returning probability of 1 neutron detected in time-window
prob_det1 = function(n1, n2, eff){
  return( prob_s1(0, n1, eff)*prob_s2(1, n2, eff) + prob_s1(1, n1, eff)*prob_s2(0, n2, eff) )
}

# Function returning probability of 2 neutron detected in time-window
prob_det2 = function(n1, n2, eff){
  return( prob_s1(2, n1, eff)*prob_s2(0, n2, eff) + 
          prob_s1(1, n1, eff)*prob_s2(1, n2, eff) +
          prob_s1(0, n1, eff)*prob_s2(2, n2, eff) )
}

# Function returning probability of 3 neutron detected in time-window
prob_det3 = function(n1, n2, eff){
  return( prob_s1(3, n1, eff)*prob_s2(0, n2, eff) + 
          prob_s1(2, n1, eff)*prob_s2(1, n2, eff) +
          prob_s1(1, n1, eff)*prob_s2(2, n2, eff) +
          prob_s1(0, n1, eff)*prob_s2(3, n2, eff)  )
}

# starting values....
# n1_mean = 193692/89424000/0.22
# n2_mean = 2625/89424000/0.22^2
n1_mean = 0.008907
dn1 =     15*0.00012
n2_mean = 623.303e-06
dn2 =     15*37e-06

norma = 13413    # Term to normalize since the numbers are VERY small!
l = 13413000/norma
m1 = 28518/norma
m2 = 459/norma
m3 = 1/norma
m0 = l-m1-m2-m3
# multiComb = factorial(l)/(factorial(m0)*factorial(m1)*factorial(m2)*factorial(m3))
tempNorm = 1/5.925e6

# Function returning observation of the whole experiment
prob_observed = function(n1, n2, eff) {
  p0 = prob_det0(n1, n2, eff)
  p1 = prob_det1(n1, n2, eff)
  p2 = prob_det2(n1, n2, eff)
  p3 = prob_det3(n1, n2, eff)
  prob_total = 1
  for (i in 1:norma) {
    prob_total = prob_total*p0^m0*p1^m1*p2^m2*p3^m3*(1-p0-p1-p2-p3)^(l-m0-m1-m2-m3)/tempNorm #*multiComb
  } 
  return( prob_total )
}

prob_observed(n1_mean, n2_mean, 0.22)

# eff vs n1
slope1 = -29.52
limit1 = 0.465
limit2 = 0.495

# n2 vs n1
slope2 = 179804.27e-6
limit4 = -750e-6
limit3 = -1000e-6

ptime <- system.time({
  
  prob_dist = data.frame(
    eff = numeric(0),     # neutron detection efficiency
    n1 = numeric(0),      # rate of room background event creating 1 neutron
    n2 = numeric(0),      # rate of room background event creating 2 neutrons
    prob = numeric(0)     # probability
  )
  for (k in 1:10000) {
    # Uniform randomly generated average rate of room background event creating 1 neutron
    n1 = runif(1, n1_mean-dn1, n1_mean+dn1) 
    prob_dist[k, "n1"] = n1
    # Uniform randomly generated average rate of room background event creating 2 neutrons
    # n2 = sample(c(n2_mean-dn2, n2_mean, n2_mean+dn2),1) 
    # n2 = runif(1, n2_mean-dn2, n2_mean+dn2)
    n2 = slope2*n1 + runif(1, limit3, limit4) 
    prob_dist[k, "n2"] = n2
    # Uniform randomly generated efficiency between 0.17% and 0.27%
    # eff = sample(c(0.217, 0.22, 0.223),1) 
    # eff = runif(1, 0.17, 0.27) 
    # eff = 0.217
    eff = slope1*n1 + runif(1, limit1, limit2) 
    prob_dist[k, "eff"] = eff
    prob_dist[k, "prob"] = prob_observed(n1, n2, eff)
  }
  
})[3]
print(ptime)

# prob_dist = prob_dist[prob_dist$prob > 1,]
# plot_ly(x = temp$n1, y = temp$eff)

eff_n1 <- plot_ly(data = prob_dist, x = ~n1, y = ~eff, color = ~prob, size = ~prob,
        text = ~paste0("Prob: ", format(prob,digits=2)),
        type = 'scatter', mode = 'markers')

n2_n1 <- plot_ly(data = prob_dist, x = ~n1, y = ~n2, color = ~prob, size = ~prob,
        text = ~paste0("Prob: ", format(prob,digits=2)),
        type = 'scatter', mode = 'markers') 

subplot( style(eff_n1, showlegend = FALSE), style(n2_n1, showlegend = FALSE), 
         nrows = 2, margin = 0.01, titleX = T, titleY = T)

plot_ly(data = prob_dist, x =~n1, y = ~eff, z = ~prob) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'n1'),
                      yaxis = list(title = 'eff'),
                      zaxis = list(title = 'prob')))

plot_ly(data = prob_dist, x = ~n1, y = ~n2, z = ~prob) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'n1'),
                      yaxis = list(title = 'n2', type = "log"),
                      zaxis = list(title = 'prob')))

fileName <- "~/Dropbox/Courses/R/alpha_n/prob_dist_room_Hendrik.RData"
save(prob_dist, file = fileName)


# Testing ------------

prob_observed(n1_mean, n2_mean, 0.2)

p1 = prob_det1(n1_mean, n2_mean, 0.2)
p1^m1

p2 = prob_det2(n1_mean, n2_mean, 0.2)
p2^m2



x <- seq(0, 5, 0.01)
plot(x, dpois(x, 2.7), ylab = "F(x)", main = "Poisson(1) CDF")

dpois(0.1,1)

# eff vs n1
x1 = c(0.00725, 0.276)
x2 = c(0.01057, 0.178)
slope1 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope1*x1[1]
slope1*0.0077+limit

# n2 vs n1
x1 = c(0.01059799, 960)
x2 = c(0.00755913, 413.6)
slope2 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope2*x1[1]
slope2*0.01+limit
