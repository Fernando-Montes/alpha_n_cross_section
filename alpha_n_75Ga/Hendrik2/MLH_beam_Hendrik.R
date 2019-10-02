# Combinatorial function
comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}

# Function returning probability of detecting det number of neutrons in time-window 
# coming from source producing orig number of neutrons
prob_det_sour = function(det, orig, rate, eff) {
  prob = 0
  for(i in ceiling(det/orig):10) {
    prob = prob + dpois(i, rate)*eff^det*(1-eff)^(orig*i-det)*comb(orig*i, det)
  }
  return( prob )
}

# Function returning probability of detected det number of neutrons in time-window
prob_det = function(det, n1, n2,  b1, b2, b3, b4, b5, effn, effb){
  probTemp = 0
  for (n1det in det:0) {
    for (n2det in (det-n1det):0) {
      for (b1det in (det-n1det-n2det):0) {
        for (b2det in (det-n1det-n2det-b1det):0) {
          for (b3det in (det-n1det-n2det-b1det-b2det):0) {
            for (b4det in (det-n1det-n2det-b1det-b2det-b3det):0) {
              b5det = det - n1det - n2det - b1det - b2det - b3det - b4det
              probTemp = probTemp + 
                prob_det_sour(n1det, 1, n1, effn)*
                prob_det_sour(n2det, 2, n2, effn)*
                prob_det_sour(b1det, 1, b1, effb)*
                prob_det_sour(b2det, 2, b2, effb)*
                prob_det_sour(b3det, 3, b3, effb)*
                prob_det_sour(b4det, 4, b4, effb)*
                prob_det_sour(b5det, 5, b5, effb)
            }
          }
        }
      }
    }
  }
  return( probTemp  )
}

# Function returning observation of the whole experiment
prob_observed = function(n1, n2, b1, b2, b3, b4, b5, effn, effb) {
  p0 = prob_det(0, n1, n2, b1, b2, b3, b4, b5, effn, effb)
  p1 = prob_det(1, n1, n2, b1, b2, b3, b4, b5, effn, effb)
  p2 = prob_det(2, n1, n2, b1, b2, b3, b4, b5, effn, effb)
  p3 = prob_det(3, n1, n2, b1, b2, b3, b4, b5, effn, effb)
  prob_total = 1
  for (i in 1:norma) {
    prob_total = prob_total*p0^m0*p1^m1*p2^m2*p3^m3*(1-p0-p1-p2-p3)^(l-m0-m1-m2-m3)/tempNorm #*multiComb
  } 
  return( prob_total )
}

norma = 8942    # Term to normalize since the numbers are VERY small!
l = 894233/norma
m1 = 2482/norma
m2 = 182/norma
m3 = 47/norma
m0 = l-m1-m2-m3  
# multiComb = factorial(l)/(factorial(m0)*factorial(m1)*factorial(m2)*factorial(m3))
tempNorm = 1/8.69
# maxSol = prob_dist[prob_dist$prob == max(prob_dist$prob),]
prob_observed(0.008113375, 0.0005280163, 0.000214365,  0, 0.0004365659, 0.0006898828, 9.642915e-05, 0.239571, 0.239571)
# prob_observed(0.007109151, 0.0003394667, 1.573719e-06, 0.00455328, 0, 0, 0, 0.2856537, 0.2856537)
# prob_observed(n1_mean, n2_mean, b1_mean, b2_mean, eff_mean)

# starting values n1, n2, eff....
fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik2/prob_dist_room_Hendrik.RData"
load(file = fileName)
prob_dist_n1_n2 = prob_dist
# plot_ly(data = prob_dist_n1_n2, x = ~n1, y = ~eff, color = ~prob, size = ~prob)
# plot_ly(data = prob_dist_n1_n2, x = ~n1, y = ~n2, color = ~prob, size = ~prob) 

# b1_mean =  4556/1271800/0.22
b1_mean = 1.2e-03
db1 = 0.4e-3
# b2_mean =  676/1271800/0.22^2
b2_mean = 0
db2 = 200e-6
b3_mean = 0
db3 = 100e-6
b4_mean = 0.75e-3
db4 = 0.20e-3
b5_mean = 0
db5 = 300e-6

# b1 vs b4
x1 = c(1076e-6, 91e-6)
x2 = c(569e-6, 255e-6)
slope1 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope1*x1[1]
slope1*x1[1]+limit
limit1_1 = 1500e-6
limit1_2 = 10000e-6

# b2 vs b4
x1 = c(406e-6, 518e-6)
x2 = c(1090e-6, 216e-6)
slope2 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope2*x1[1]
slope2*x1[1]+limit
limit2_1 = 400e-6
limit2_2 = 700e-6

# b3 vs b4
x1 = c(406e-6, 241e-6)
x2 = c(1075e-6, 62e-6)
slope3 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope3*x1[1]
slope3*x1[1]+limit
limit3_1 = 0e-6
limit3_2 = 500e-6

# b4 vs effb
x1 = c(0.277, 687e-6)
x2 = c(0.246, 761e-6)
slope4 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope4*x1[1]
slope4*x1[1]+limit
limit4_1 = 0.0013
limit4_2 = 0.0025

ptime <- system.time({
  
  prob_dist = data.frame(
    effn = numeric(0),    # neutron detection efficiency room background neutrons
    effb = numeric(0),    # neutron detection efficiency beam-induced background neutrons
    n1 = numeric(0),      # rate of room background event creating 1 neutron
    n2 = numeric(0),      # rate of room background event creating 2 neutrons
    b1 = numeric(0),      # rate of beam-induced background event creating 1 neutron
    b2 = numeric(0),      # rate of beam-induced background event creating 2 neutrons
    b3 = numeric(0),      # rate of beam-induced background event creating 3 neutrons
    b4 = numeric(0),      # rate of beam-induced background event creating 4 neutrons
    b5 = numeric(0),      # rate of beam-induced background event creating 5 neutrons
    prob = numeric(0)     # probability
  )
  ind = rownames(prob_dist_n1_n2[prob_dist_n1_n2$eff>0.17,])
  for (k in 1:10000) {
    # Quantities obtained from MLH analysis of n1 and n2 data
    indx = sample(ind,1)
    effn = prob_dist_n1_n2[indx,"eff"]
    n1 =   prob_dist_n1_n2[indx,"n1"]
    n2 =   prob_dist_n1_n2[indx,"n2"]
    prob_dist[k, "n1"] = n1
    prob_dist[k, "n2"] = n2
    prob_dist[k, "effn"] = effn
    
    # New quantities to be determined:
    # effb = slope4*b4 + runif(1, limit4_1, limit4_2) 
    # if (effb<0) effb=0
    effb = runif(1, effn - 0.05, effn + 0.10)
    # effb = effn
    prob_dist[k, "effb"] = effb
    # b4 = slope4*effb + runif(1, limit4_1, limit4_2) 
    # if (b4<0) b4=0
    b4 = runif(1, 0.0005, 0.0009)
    prob_dist[k, "b4"] = b4
    # b1 = slope1*b4 + runif(1, limit1_1, limit1_2)
    # if (b1<0) b1=0
    b1 = runif(1, b1_mean-db1, b1_mean+db1)
    # b1 = b1_mean*10^runif(1, -2, 0)
    # b1 = 0
    prob_dist[k, "b1"] = b1
    # b2 = slope2*b4 + runif(1, limit2_1, limit2_2)
    # if (b2<0) b2=0
    # b2 = runif(1, 0, b2_mean+db2)
    b2 = 0
    prob_dist[k, "b2"] = b2
    # b3 = slope3*b4 + runif(1, limit3_1, limit3_2) 
    # if (b3<0) b3=0
    b3 = runif(1, 0, b3_mean+db3)
    # b3 = 0 
    prob_dist[k, "b3"] = b3
    b5 = runif(1, 0, b5_mean+db5)
    # b5 = 0
    prob_dist[k, "b5"] = b5
    prob_dist[k, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, effn, effb)
  }
  
})[3]
print(ptime)

b1b4 <- plot_ly(data = prob_dist, x = ~b4, y = ~b1, color = ~prob, size = ~prob, 
        text = ~paste0("Prob: ", format(prob,digits=2)),
        type = 'scatter', mode = 'markers')
b2b4 <- plot_ly(data = prob_dist, x = ~b4, y = ~b2, color = ~prob, size = ~prob, 
                text = ~paste0("Prob: ", format(prob,digits=2)),
                type = 'scatter', mode = 'markers') 
b3b4 <- plot_ly(data = prob_dist, x = ~b4, y = ~b3, color = ~prob, size = ~prob, 
                text = ~paste0("Prob: ", format(prob,digits=2)),
                type = 'scatter', mode = 'markers') 
b5b4 <- plot_ly(data = prob_dist, x = ~b4, y = ~b5, color = ~prob, size = ~prob, 
                text = ~paste0("Prob: ", format(prob,digits=2)),
                type = 'scatter', mode = 'markers') 
b4effb <- plot_ly(data = prob_dist, x = ~effb, y = ~b4, color = ~prob, size = ~prob, 
                text = ~paste0("Prob: ", format(prob,digits=2)),
                type = 'scatter', mode = 'markers') 
effbeffn <- plot_ly(data = prob_dist, x = ~effn, y = ~effb, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers') 

sb1 <- subplot( style(b1b4, showlegend = FALSE), style(b2b4, showlegend = FALSE), style(b3b4, showlegend = FALSE),
        nrows = 3, margin = 0.01, titleX = T, titleY = T)
sb2 <- subplot( style(b5b4, showlegend = FALSE), style(b4effb, showlegend = FALSE), style(effbeffn, showlegend = FALSE), 
        nrows = 3, margin = 0.01, titleX = T, titleY = T)

subplot(sb1, sb2, margin = 0.05, titleX = T, titleY = T)

# 3D plot ----------------
plot_ly(data = prob_dist, x = ~b2, y = ~b3, z = ~prob, color = ~prob, size = ~prob ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'b2', type = "log"),
                      yaxis = list(title = 'b3'),
                      zaxis = list(title = 'prob')))

fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik2/prob_dist_beam_Hendrik.RData"
save(prob_dist, file = fileName)

# Testing ------------
maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
prob_dist = prob_dist[prob_dist$prob > maxP$prob*0.1,]

keep2 = prob_dist[prob_dist$prob == max(prob_dist$prob),]

prob_observed(n1_mean, n2_mean, 0.2)

p0 = prob_det(0, n1_mean, n2_mean, b1_mean, b2_mean, eff_mean)
p1 = prob_det(1, n1_mean, n2_mean, b1_mean, b2_mean, eff_mean)

prob_total = 1

  prob_total = prob_total*p0^m0*p1^m1*p2^m2*p3^m3*(1-p0-p1-p2-p3)^(l-m0-m1-m2-m3)/tempNorm #*multiComb


x <- seq(0, 5, 0.01)
plot(x, dpois(x, 2.7), ylab = "F(x)", main = "Poisson(1) CDF")

dpois(0.1,1)

# b1 vs b4
x1 = c(568e-6, 339e-6)
x2 = c(3400e-6, 0)
slope1 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope1*x1[1]
slope1*x1[1]+limit
slope1*1.893487e-03+limit
limit1_1 = 0
limit1_2 = 800e-6

# b2 vs b4
x1 = c(568e-6, 654e-6)
x2 = c(3400e-6, 0)
slope2 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope2*x1[1]
slope2*x1[1]+limit
limit2_1 = 0
limit2_2 = 1500e-6

# b3 vs b4
x1 = c(568e-6, 1721e-6)
x2 = c(3500e-6, 0)
slope3 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope3*x1[1]
slope3*x1[1]+limit
limit3_1 = 1500e-6
limit3_2 = 2500e-6

# effb vs b4
x1 = c(568e-6, 0.286)
x2 = c(3500e-6, 0.185)
slope4 = (x2[2]-x1[2])/(x2[1]-x1[1])
limit = x1[2]-slope4*x1[1]
slope4*x1[1]+limit
limit4_1 = 0.270
limit4_2 = 0.330
