# Parallel computing
library(doParallel)
registerDoParallel(cores=4)
library(plotly)

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

norma = 2000    # Term to normalize since the numbers are VERY small!
l = 568830/norma
m1 = 2147/norma
m2 = 315/norma
m3 = 31/norma
m0 = l-m1-m2-m3
# multiComb = factorial(l)/(factorial(m0)*factorial(m1)*factorial(m2)*factorial(m3))
tempNorm = 1/5e3
prob_observed(4.571312e-05, 1.205991e-06, 0.003700, 0.000550, 0.000050, 0, 0, 1.0, 1.0)

# starting values n1, n2, eff....
fileName <- "~/Dropbox/Courses/R/alpha_n/alpha_n_85Br/3.94MeV/prob_dist_room.RData"
load(file = fileName)
prob_dist_n1_n2 = prob_dist
# plot_ly(data = prob_dist_n1_n2, x = ~n1, y = ~eff, color = ~prob, size = ~prob)
# plot_ly(data = prob_dist_n1_n2, x = ~n1, y = ~n2, color = ~prob, size = ~prob) 

# b1_mean =  4556/1271800/0.22
b1_mean = 3750e-6
db1 = 200e-06
b2_mean = 550e-6
db2 = 100e-6
b3_mean = 50e-6
db3 = 30e-6
# b4_mean = 0.75e-3
# db4 = 0.20e-3
# b5_mean = 0
# db5 = 300e-6

ptime <- system.time({
  ind = rownames(prob_dist_n1_n2)
  table <- foreach(k=1:1000, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
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
    # Quantities obtained from MLH analysis of n1 and n2 data
    indx = sample(ind,1)
    effn = prob_dist_n1_n2[indx,"eff"]
    n1 =   prob_dist_n1_n2[indx,"n1"]
    n2 =   prob_dist_n1_n2[indx,"n2"]
    prob_dist[1, "n1"] = n1
    prob_dist[1, "n2"] = n2
    prob_dist[1, "effn"] = effn
    
    # New quantities to be determined:
    # effb = runif(1, effn - 0.05, effn + 0.10)
    effb = 1.0
    prob_dist[1, "effb"] = effb
    b1 = runif(1, b1_mean-db1, b1_mean+db1)
    # b1 = runif(1, 0, db1)
    # b1 = 0
    prob_dist[1, "b1"] = b1
    b2 = runif(1, b2_mean-db2, b2_mean+db2)
    # b2 = 0
    prob_dist[1, "b2"] = b2
    # b3 = slope3*b4 + runif(1, limit3_1, limit3_2) 
    # if (b3<0) b3=0
    b3 = runif(1, b3_mean-db3, b3_mean+db3)
    # b3 = 0 
    prob_dist[1, "b3"] = b3
    # b4 = runif(1, 0.0005, 0.0009)
    b4 = 0
    prob_dist[1, "b4"] = b4
    # b5 = runif(1, 0, b5_mean+db5)
    b5 = 0
    prob_dist[1, "b5"] = b5
    prob_dist[1, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, effn, effb)
    # if (prob_dist[1, "prob"] > 1e+200) return ( prob_dist )
    return ( prob_dist )
  }
  prob_dist = table
  
})[3]
print(ptime)

b1b2 <- plot_ly(data = prob_dist, x = ~b1, y = ~b2, color = ~prob, size = ~prob, 
        text = ~paste0("Prob: ", format(prob,digits=2)),
        type = 'scatter', mode = 'markers')

b1b3 <- plot_ly(data = prob_dist, x = ~b1, y = ~b3, color = ~prob, size = ~prob, 
        text = ~paste0("Prob: ", format(prob,digits=2)),
        type = 'scatter', mode = 'markers')

subplot(b1b2, b1b3, margin = 0.05, nrows=2,  titleX = T, titleY = T)


# 3D plot ----------------
plot_ly(data = prob_dist, x = ~b2, y = ~b3, z = ~prob, color = ~prob, size = ~prob ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'b2', type = "log"),
                      yaxis = list(title = 'b3'),
                      zaxis = list(title = 'prob')))

fileName <- "~/Dropbox/Courses/R/alpha_n/alpha_n_85Br/3.94MeV/prob_dist_beam.RData"
save(prob_dist, file = fileName)

# Testing ------------
maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
prob_dist = prob_dist[prob_dist$prob > maxP$prob*0.5,]

