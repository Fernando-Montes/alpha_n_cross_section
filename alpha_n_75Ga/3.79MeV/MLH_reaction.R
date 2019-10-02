library(plotly)
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

beam_ratio = 1.087463 # ratio of beam intensity between gas in and gas out runs
norma = 9854    # Term to normalize since the numbers are VERY small!
l = 1256200/norma
m1 = 6390/norma
m2 = 874/norma
m3 = 312/norma
m0 = l-m1-m2-m3
# multiComb = factorial(l)/(factorial(m0)*factorial(m1)*factorial(m2)*factorial(m3))
factorCS = 1256200/(1768844*51.47)/(1.36022E+25*0.0399796)*1E+31 # convert rate to mbarn
tempNorm = 1/1.655e2
prob_observed(0.002272948, 2.955637e-05, 0.001838473, 0.0005903588, 9.590115e-05, 0, 0, 0.004133199, 0.0001940544, 1.0, 1.0, 0.2396504)

# starting values n1, n2, eff....
fileName <- "~/Dropbox/Courses/R/alpha_n/3.79MeV/prob_dist_beam.RData"
load(file = fileName)
prob_dist_beam = prob_dist
# plot_ly(data = prob_dist_beam, x = ~n1, y = ~eff, color = ~prob, size = ~prob)
# plot_ly(data = prob_dist_beam, x = ~n1, y = ~n2, color = ~prob, size = ~prob) 

a1_mean = 1150/factorCS
da1 =     700/factorCS
a2_mean = 200/factorCS
da2 =     200/factorCS

ptime <- system.time({
  
  ind = rownames(prob_dist_beam)
  table <- foreach(k=1:1000000, .combine = rbind, .verbose = F, .errorhandling = "remove") %dopar% { 
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
    a1 = runif(1, a1_mean-da1, a1_mean+da1)
    # a1 = 0
    prob_dist[1, "a1"] = a1
    a2 = runif(1, a2_mean-da2, a2_mean+da2)
    # a2 = 0
    prob_dist[1, "a2"] = a2
    effr = runif(1, 0.17, 0.27)
    #effr = 0.22
    prob_dist[1, "effr"] = effr
    
    prob_dist[1, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
    if (prob_dist[1, "prob"] > 1e0) return ( prob_dist )
    # return ( prob_dist )
  }
  
})[3]
print(ptime)
prob_dist = table


prob_dist = rbind(prob_dist, prob_dist_cum)
prob_dist_cum = prob_dist


a2a1CS <- plot_ly(data = prob_dist, x = ~a1*factorCS, y = ~a2*factorCS, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers')

a1effr <- plot_ly(data = prob_dist, x = ~effr, y = ~a1*factorCS, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers')
a2effr <- plot_ly(data = prob_dist, x = ~effr, y = ~a2*factorCS, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers')

subplot( style(a2a1CS, showlegend = FALSE), style(a1effr, showlegend = FALSE), style(a2effr, showlegend = FALSE),
         nrows = 3, margin = 0.01, titleX = T, titleY = T)



# 3D plot ----------------
plot_ly(data = prob_dist, x = ~a1*factorCS, y = ~a2*factorCS, z = ~prob, color = ~prob, size = ~prob, 
        type="scatter3d", mode='markers' ) %>%
  layout(scene = list(xaxis = list(title = 'a1'),
                      yaxis = list(title = 'a2'),
                      zaxis = list(title = 'prob')))

fileName <- "~/Dropbox/Courses/R/alpha_n/3.79MeV/prob_dist_reaction.RData"
save(prob_dist, file = fileName)

# Testing ------------
maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
prob_dist = prob_dist[prob_dist$prob > maxP$prob*0.001,]




# Confidence intervals -----------
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


plot_ly(z = prob_table, y=rownames(prob_table), x=colnames(prob_table), type = "heatmap") 

plot_ly(z = prob_table, y=rownames(prob_table), x=colnames(prob_table), type = "surface")

plot_ly(
  z = prob_table,
  y = rownames(prob_table),
  x = colnames(prob_table),
  type = "contour",
  contours = list(showlabels = FALSE, coloring = 'heatmap')
)  %>%
  colorbar(title = "Probability") %>%
  layout(xaxis = list(title = "a1n"),
         yaxis = list(title = "a2n"))


z = t(prob_table)
y = as.numeric ( rownames(prob_table) ) #a2n
x = as.numeric ( colnames(prob_table) ) #a1n

uniqueVals <- rev(unique(sort(z)))
cumProbs <- sapply(uniqueVals,
                   function(val) sum( z[z>=val] ) )
names(cumProbs) <- uniqueVals
cumDensity <- matrix(nrow=nrow(z), ncol=ncol(z))
cumDensity[] <- cumProbs[as.character(z)]

library(grDevices)
plot(0:2000, 0:2000, type = "n", xlim=c(600, 2000), ylim=c(0, 300), xlab="a1n [mb]", ylab="a2n [mb]", axes = F)  # setting up coord. system
axis(side = 1, at = c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000))
axis(side = 2, at = c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000))
box()
points(x = maxP$a1*factorCS, y = maxP$a2*factorCS, col = "blue")
# points(x = 34, y = 780, col = "blue")
contour(x, y, cumDensity,
        xlab="a1n [mb]", ylab="a2n [mb]",
        levels=c(0.59, .80), lwd=4, lty=2, add = T, col = "pink")

# ---------------------------------------------------
# Confidence intervals for the sum (a,xn) -----------
# Create cuts ---
axn <- cut((prob_dist$a1*factorCS+prob_dist$a2*factorCS), 30)
prob_table_axn <- table(axn)
# resetting probability table
for (i in 1:dim(prob_table_axn)) {
    prob_table_axn[i] = 0
}
# filling probabilities
for (i in 1:dim(prob_dist)[1]) {
  prob_table_axn[axn[i]] = prob_table_axn[axn[i]] + prob_dist$prob[i]
}
# resetting rownames
for (i in 1:length(rownames(prob_table_axn))) {
  t1 = unlist( strsplit(rownames(prob_table_axn)[i],",") )
  t2 = as.numeric(unlist( strsplit(t1[1],"[(]") )[2])
  t3 = as.numeric(unlist( strsplit(t1[2],"]") )[1])
  rownames(prob_table_axn)[i] = as.character((t2+t3)/2)
}
# integrating
totalsum = 0
for (i in 1:dim(prob_table_axn)) {
    # if (prob_table[i,j] > totalsum) totalsum = prob_table[i,j]
    totalsum = totalsum + prob_table_axn[i]
}
# normalizing
for (i in 1:dim(prob_table_axn)) {
  prob_table_axn[i] = prob_table_axn[i]/totalsum
}

plot_ly(y = prob_table_axn, x=rownames(prob_table_axn), type = "bar") 

max(prob_table_axn)
prob_table_axn[14]
sum(prob_table_axn[1:9])
prob_table_axn[9]
sum(prob_table_axn[20:30])
prob_table_axn[20]


# # First present the data in a data-frame
# tab <- as.data.frame(prob_table_axn)
# #Apply function nls
# (res <- nls( Freq ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=1500,sigma=100,k=0.1) , data = tab))
# 
# v <- summary(res)$parameters[,"Estimate"]
# plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2), xlim=c(790,2000) , col = "blue")
# plot(x = tab$axn, y=tab$Freq, col = "pink")