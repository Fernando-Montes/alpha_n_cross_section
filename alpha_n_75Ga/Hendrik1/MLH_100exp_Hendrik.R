library(plotly)
library(gplots)
# Parallel computing
library(doParallel)
registerDoParallel(cores=4)

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
prob_det = function(det, n1, n2,  b1, b2, b3, b4, b5, a1, a2, effn, effb, effr){
  probTemp = 0
  for (n1det in det:0) {
    for (n2det in (det-n1det):0) {
      for (b1det in (det-n1det-n2det):0) {
        for (b2det in (det-n1det-n2det-b1det):0) {
          for (b3det in (det-n1det-n2det-b1det-b2det):0) {
            for (b4det in (det-n1det-n2det-b1det-b2det-b3det):0) {
              for (b5det in (det-n1det-n2det-b1det-b2det-b3det-b4det):0) {
                for (a1det in (det-n1det-n2det-b1det-b2det-b3det-b4det-b5det):0) {
                  a2det = det - n1det - n2det - b1det - b2det - b3det - b4det - b5det - a1det
                  probTemp = probTemp + 
                    prob_det_sour(n1det, 1, n1, effn)*
                    prob_det_sour(n2det, 2, n2, effn)*
                    prob_det_sour(b1det, 1, b1, effb)*
                    prob_det_sour(b2det, 2, b2, effb)*
                    prob_det_sour(b3det, 3, b3, effb)*
                    prob_det_sour(b4det, 4, b4, effb)*
                    prob_det_sour(b5det, 5, b5, effb)*
                    prob_det_sour(a1det, 1, a1, effr)*
                    prob_det_sour(a2det, 2, a2, effr)
                }
              }
            }
          }
        }
      }
    }
  }
  return( probTemp  )
}

# Function returning observation of the whole experiment
prob_observed = function(n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr) {
  p0 = prob_det(0, n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
  p1 = prob_det(1, n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
  p2 = prob_det(2, n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
  p3 = prob_det(3, n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
  prob_total = 1
  for (i in 1:norma) {
    prob_total = prob_total*p0^m0*p1^m1*p2^m2*p3^m3*(1-p0-p1-p2-p3)^(l-m0-m1-m2-m3)/tempNorm #*multiComb
  } 
  return( prob_total )
}

fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik/fernando_all_gasin.csv"
gas_in = read.csv(file = fileName, skip = 2, header = FALSE)
fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik/fernando_all_gasout.csv"
gas_out = read.csv(file = fileName, skip = 2, header = FALSE)

for (jj in 21:21) {
  print(jj)
  rm(prob_dist_n1_n2)
  rm(prob_dist_beam)
  rm(prob_dist)
  
  print("starting room background")
  # ---------------------
  # Room background -----
  # ---------------------
  norma = 13413    # Term to normalize since the numbers are VERY small!
  l = (13413000-gas_in[jj,5])/norma
  m1 = gas_in[jj,2]/norma
  m2 = gas_in[jj,3]/norma
  m3 = gas_in[jj,4]/norma
  m0 = l-m1-m2-m3
  
  norma = 13413    # Term to normalize since the numbers are VERY small!
  l = (13413500-0)/norma
  m1 = 28650/norma
  m2 = 490/norma
  m3 = 2/norma
  m0 = l-m1-m2-m3
  
  # region ----
  n1_mean = 0.008907
  dn1 =     15*0.00012
  n2_mean = 623.303e-06
  dn2 =     15*37e-06
  # tempNorm = 1/5.655e6
  tempNorm = 1/6.46e6
  if ( prob_observed(0.007209895, 0.0004334408, 0, 0, 0, 0, 0, 0, 0, 0.2708555, 0, 0) < 1e-300 ) { break }
  if ( prob_observed(0.007209895, 0.0004334408, 0, 0, 0, 0, 0, 0, 0, 0.2708555, 0, 0) > 1e+200 ) { break }
  
  # eff vs n1 ----
  slope1 = -29.52
  limit1 = 0.465
  limit2 = 0.495
  
  # n2 vs n1 ----
  slope2 = 179804.27e-6
  limit4 = -750e-6
  limit3 = -1000e-6
  
  # loop ----
  table <- foreach(k=1:10000, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
    prob_dist = data.frame(
      eff = numeric(0),     # neutron detection efficiency
      n1 = numeric(0),      # rate of room background event creating 1 neutron
      n2 = numeric(0),      # rate of room background event creating 2 neutrons
      prob = numeric(0)     # probability
    )
    # Uniform randomly generated average rate of room background event creating 1 neutron
    n1 = runif(1, n1_mean-dn1, n1_mean+dn1) 
    prob_dist[1, "n1"] = n1
    # Uniform randomly generated average rate of room background event creating 2 neutrons
    # n2 = sample(c(n2_mean-dn2, n2_mean, n2_mean+dn2),1) 
    # n2 = runif(1, n2_mean-dn2, n2_mean+dn2)
    n2 = slope2*n1 + runif(1, limit3, limit4) 
    prob_dist[1, "n2"] = n2
    # Uniform randomly generated efficiency between 0.17% and 0.27%
    # eff = sample(c(0.217, 0.22, 0.223),1) 
    # eff = runif(1, 0.17, 0.27) 
    # eff = 0.217
    eff = slope1*n1 + runif(1, limit1, limit2) 
    prob_dist[1, "eff"] = eff
    prob_dist[1, "prob"] = prob_observed(n1, n2, 0, 0, 0, 0, 0, 0, 0, eff, 0, 0)
    return ( prob_dist )
  }
  prob_dist = table
  
  maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
  prob_dist_n1_n2 = prob_dist[prob_dist$prob > maxP$prob*0.1,]
  
  print("starting beam background")
  # ---------------------
  # Beam background -----
  # ---------------------
  norma = 8942    # Term to normalize since the numbers are VERY small!
  l = (894200-gas_out[jj,10])/norma
  m1 = gas_out[jj,7]/norma
  m2 = gas_out[jj,8]/norma
  m3 = gas_out[jj,9]/norma
  m0 = l-m1-m2-m3  
  
  norma = 9121    # Term to normalize since the numbers are VERY small!
  l = (912118-4)/norma
  m1 = 2372/norma
  m2 = 195/norma
  m3 = 51/norma
  m0 = l-m1-m2-m3  
  
  # region ----
  b1_mean = 1.2e-03
  db1 = 0.4e-3
  b2_mean = 0
  db2 = 200e-6
  b3_mean = 0
  db3 = 100e-6
  b4_mean = 0.75e-3
  db4 = 0.20e-3
  b5_mean = 0
  db5 = 300e-6
  tempNorm = 1/7.93
  if (prob_observed(0.008312573, 0.0005761838, 0.0009612136,  0, 3.541429e-05, 0.0007501725, 2.341676e-06, 0, 0, 0.2325052, 0.2677564, 0) < 1e-200) { break }
  if (prob_observed(0.008312573, 0.0005761838, 0.0009612136,  0, 3.541429e-05, 0.0007501725, 2.341676e-06, 0, 0, 0.2325052, 0.2677564, 0) > 1e+200) {break}
    
  # b1 vs b4 ----
  x1 = c(1076e-6, 91e-6)
  x2 = c(569e-6, 255e-6)
  slope1 = (x2[2]-x1[2])/(x2[1]-x1[1])
  limit = x1[2]-slope1*x1[1]
  slope1*x1[1]+limit
  limit1_1 = 1500e-6
  limit1_2 = 10000e-6
  
  # b2 vs b4 ----
  x1 = c(406e-6, 518e-6)
  x2 = c(1090e-6, 216e-6)
  slope2 = (x2[2]-x1[2])/(x2[1]-x1[1])
  limit = x1[2]-slope2*x1[1]
  slope2*x1[1]+limit
  limit2_1 = 400e-6
  limit2_2 = 700e-6
  
  # b3 vs b4 ----
  x1 = c(406e-6, 241e-6)
  x2 = c(1075e-6, 62e-6)
  slope3 = (x2[2]-x1[2])/(x2[1]-x1[1])
  limit = x1[2]-slope3*x1[1]
  slope3*x1[1]+limit
  limit3_1 = 0e-6
  limit3_2 = 500e-6
  
  # b4 vs effb ----
  x1 = c(0.277, 687e-6)
  x2 = c(0.246, 761e-6)
  slope4 = (x2[2]-x1[2])/(x2[1]-x1[1])
  limit = x1[2]-slope4*x1[1]
  slope4*x1[1]+limit
  limit4_1 = 0.0013
  limit4_2 = 0.0025
  
  # loop ----
  ind = rownames(prob_dist_n1_n2)
  table <- foreach(k=1:10000, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
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
    effb = runif(1, effn - 0.05, effn + 0.05)
    prob_dist[1, "effb"] = effb
    b4 = runif(1, 0.0005, 0.0009)
    prob_dist[1, "b4"] = b4
    b1 = runif(1, b1_mean-db1, b1_mean+db1)
    prob_dist[1, "b1"] = b1
    b2 = 0
    prob_dist[1, "b2"] = b2
    b3 = runif(1, 0, b3_mean+db3)
    prob_dist[1, "b3"] = b3
    b5 = runif(1, 0, b5_mean+db5)
    prob_dist[1, "b5"] = b5
    prob_dist[1, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, 0, 0, effn, effb, 0)
    return ( prob_dist )
  }
  prob_dist = table
  
  maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
  prob_dist_beam = prob_dist[prob_dist$prob > maxP$prob*0.1,]  

  

  print("starting reaction analysis")
  # ---------------------
  # Reaction ------------
  # ---------------------
  beam_ratio = 1 # ratio of beam intensity between gas in and gas out runs
  norma = 8942    # Term to normalize since the numbers are VERY small!
  l = (894200-gas_in[jj,10])/norma
  m1 = gas_in[jj,7]/norma
  m2 = gas_in[jj,8]/norma
  m3 = gas_in[jj,9]/norma
  m0 = l-m1-m2-m3
  
  norma = 9121    # Term to normalize since the numbers are VERY small!
  l = (912118-2)/norma
  m1 = 4169/norma
  m2 = 210/norma
  m3 = 44/norma
  m0 = l-m1-m2-m3
  
  # region ----
  factorCS = 1/(5.43e19*2000/0.005*0.0003)*10^27
  a1_mean = 6.0e-03
  da1 =     1.5e-3
  a2_mean = 0.9e-3
  da2 =     0.4e-3
  tempNorm = 1/24.23
  if (prob_observed(0.008681848, 0.0007058741, 0.001282992,  0, 1.449724e-05, 0.0007181153, 2.26838e-05, 0.009353153, 0.0001608321,
                    0.2174893, 0.2647424, 0.1837667) < 1e-160) { break }
  if (prob_observed(0.008681848, 0.0007058741, 0.001282992,  0, 1.449724e-05, 0.0007181153, 2.26838e-05, 0.009353153, 0.0001608321,
                         0.2174893, 0.2647424, 0.1837667) > 1e+200) { break }
  
  
  # loop ----
  ind = rownames(prob_dist_beam)
  table <- foreach(k=1:500000, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
    prob_dist = data.frame(
      effn = numeric(0),    # neutron detection efficiency room background neutrons
      effb = numeric(0),    # neutron detection efficiency beam-induced background neutrons
      effr = numeric(0),    # neutron detection efficiency reaction-induced neutrons
      n1 = numeric(0),      # rate of room background event creating 1 neutron
      n2 = numeric(0),      # rate of room background event creating 2 neutrons
      b1 = numeric(0),      # rate of beam-induced background event creating 1 neutron
      b2 = numeric(0),      # rate of beam-induced background event creating 2 neutrons
      b3 = numeric(0),      # rate of beam-induced background event creating 3 neutrons
      b4 = numeric(0),      # rate of beam-induced background event creating 4 neutrons
      b5 = numeric(0),      # rate of beam-induced background event creating 5 neutrons
      a1 = numeric(0),      # rate of reaction-induced event creating 1 neutron
      a2 = numeric(0),      # rate of reaction-induced event creating 2 neutrons
      prob = numeric(0)     # probability
    )
    # Quantities obtained from MLH analysis of room and beam-induced background data
    indx = sample(ind,1)
    print(indx)
    effn = prob_dist_beam[indx,"effn"]
    effb = prob_dist_beam[indx,"effb"]
    effr = prob_dist_beam[indx,"effr"]
    n1 =   prob_dist_beam[indx,"n1"]
    n2 =   prob_dist_beam[indx,"n2"]
    b1 =   prob_dist_beam[indx,"b1"]*beam_ratio  
    b2 =   prob_dist_beam[indx,"b2"]*beam_ratio
    b3 =   prob_dist_beam[indx,"b3"]*beam_ratio
    b4 =   prob_dist_beam[indx,"b4"]*beam_ratio
    b5 =   prob_dist_beam[indx,"b5"]*beam_ratio
    prob_dist[1, "effn"] = effn
    prob_dist[1, "effb"] = effb
    prob_dist[1, "n1"] = n1
    prob_dist[1, "n2"] = n2
    prob_dist[1, "b1"] = b1
    prob_dist[1, "b2"] = b2
    prob_dist[1, "b3"] = b3
    prob_dist[1, "b4"] = b4
    prob_dist[1, "b5"] = b5
    
    # New quantities to be determined:
    a1 = runif(1, 0.005, 0.012)
    # a1 = 0
    prob_dist[1, "a1"] = a1
    a2 = runif(1, 0, a2_mean+da2)
    # a2 = 0
    prob_dist[1, "a2"] = a2
    effr = runif(1, 0.17, 0.27)
    prob_dist[1, "effr"] = effr
    
    prob_dist[1, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
    return ( prob_dist )
  }
  prob_dist = table
  
  maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
  prob_dist = prob_dist[prob_dist$prob > maxP$prob*0.001,]
  
  # Confidence intervals ----
   
  # Create cuts ---
  a1_c <- cut(prob_dist$a1*factorCS, 30)
  a2_c <- cut(prob_dist$a2*factorCS, 30)
  prob_table <- table(a2_c, a1_c)
  # resetting probability table
  for (i in 1:dim(prob_table)[1]) {
    for (j in 1:dim(prob_table)[2]) {
      prob_table[i,j] = 0
    }
  }
  for (i in 1:dim(prob_dist)[1]) {
    prob_table[a2_c[i], a1_c[i]] = prob_table[a2_c[i], a1_c[i]] + prob_dist$prob[i]
  }
  # resetting colnames and rownames
  for (i in 1:length(colnames(prob_table))) {
    t1 = unlist( strsplit(colnames(prob_table)[i],",") )
    t2 = as.numeric(unlist( strsplit(t1[1],"[(]") )[2])
    t3 = as.numeric(unlist( strsplit(t1[2],"]") )[1])
    colnames(prob_table)[i] = as.character((t2+t3)/2)
  }
  for (i in 1:length(rownames(prob_table))) {
    t1 = unlist( strsplit(rownames(prob_table)[i],",") )
    t2 = as.numeric(unlist( strsplit(t1[1],"[(]") )[2])
    t3 = as.numeric(unlist( strsplit(t1[2],"]") )[1])
    rownames(prob_table)[i] = as.character((t2+t3)/2)
  }
  # integrating
  totalsum = 0
  for (i in 1:dim(prob_table)[1]) {
    for (j in 1:dim(prob_table)[2]) {
      # if (prob_table[i,j] > totalsum) totalsum = prob_table[i,j]
      totalsum = totalsum + prob_table[i,j]
    }
  }
  # normalizing
  for (i in 1:dim(prob_table)[1]) {
    for (j in 1:dim(prob_table)[2]) {
      prob_table[i,j] = prob_table[i,j]/totalsum
    }
  }

  # plot_ly(z = prob_table, y=rownames(prob_table), x=colnames(prob_table), type = "contour",
  #         contours = list( showlabels = TRUE, coloring = 'heatmap' )
  #         )  %>%
  #   colorbar(title = "Probability") %>%
  #   layout(xaxis = list(title = "a1n"), yaxis = list(title = "a2n"))

  
  # Using Ci2d -----
  # let ci2d handle plotting contours...
  # est <- ci2d(prob_dist$a1*factorCS, prob_dist$a2*factorCS, show="contour", col="black", show.points = F, 
  #             xlab="a1n [mb]", ylab="a2n [mb]",
  #             ci.levels=c(0.39, 0.68, 0.90, 0.95))
  # est <- ci2d(prob_dist$a1*factorCS, prob_dist$a2*factorCS, show="contour", col="black", show.points = F, 
  #             xlab="a1n [mb]", ylab="a2n [mb]",
  #             ci.levels=c(0.39, 0.68, 0.90, 0.95))
  
  z = t(prob_table)
  y = as.numeric ( rownames(prob_table) ) #a2n
  x = as.numeric ( colnames(prob_table) ) #a1n
  
  uniqueVals <- rev(unique(sort(z)))
  cumProbs <- sapply(uniqueVals,
                     function(val) sum( z[z>=val] ) )
  names(cumProbs) <- uniqueVals
  cumDensity <- matrix(nrow=nrow(z), ncol=ncol(z))
  cumDensity[] <- cumProbs[as.character(z)]

  contour(x, y, cumDensity,
          xlab="a1n [mb]", ylab="a2n [mb]",
          levels=c(0.59, .80), lwd=4, lty=2)
  
  aa1n = 1459
  aa2n = 55.7
  res = cumDensity[which(abs((x-aa1n)) == min(abs(x-aa1n))), which(abs((y-aa2n)) == min(abs(y-aa2n)))]
  
  # fileNameRes <- "~/Dropbox/Courses/R/alpha_n/Hendrik/ResultErrorBar.csv"
  # write(res, file=fileNameRes, append=TRUE)
}

# Estimating contour -----
fileNameRes <- "~/Dropbox/Courses/R/alpha_n/Hendrik/ResultErrorBar.csv"
temp = read.csv(file = fileNameRes, skip = 0, header = FALSE)

res = data.frame(
  vol = numeric(0),    
  conf = numeric(0)
)

j = 1
for ( i in seq(0,1,0.05) ) {
  res[j,1] = i
  res[j,2] = sum(i > temp)/dim(temp)[1]
  j = j+1
}

plot_ly(data = res, x = ~vol, y = ~conf, 
        type = 'scatter', mode = 'markers') %>% 
  layout(  xaxis = list(title = "Contour [%]"), yaxis = list(title = "Confidence [%]") )

# Fitting function with a neural network
# Load NN package
library(RSNNS)
# Fitting with a Neural Network
model <- mlp(res$vol, res$conf, size=c(3,3,3), 
             maxit=100000, linOut=TRUE, 
             learnFuncParams=c(0.001, 0.001), 
             hiddenActFunc="Act_TanH")

y_pred <- predict(model, as.matrix(seq(0,1,0.01)))
resfit = as.data.frame( cbind(seq(0,1,0.01), y_pred) )
colnames(resfit) = c("vol", "fit")

plot_ly(data = res, x = ~vol, y = ~conf, 
        type = 'scatter', mode = 'markers')  %>% 
  add_lines(data = resfit, x = ~vol, y = ~fit, 
            line = list(color = 'black'), showlegend = FALSE) %>% 
  layout(  xaxis = list(title = "Contour [%]"), yaxis = list(title = "Confidence [%]") )

predict(model, 0.67)