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

beam_ratio = 1 # ratio of beam intensity between gas in and gas out runs
norma = 88529    # Term to normalize since the numbers are VERY small!
l = (88529067-266)/norma
m1 = 384041/norma
m2 = 36059/norma
m3 = 5010/norma
m0 = l-m1-m2-m3
# multiComb = factorial(l)/(factorial(m0)*factorial(m1)*factorial(m2)*factorial(m3))
factorCS = 985400/(1271948*51.47)/5.43159E+23*1E+31 # convert rate to mbarn (experiment)
factorCS = 1/(5.43e19*2000/0.005*0.0003)*10^27 # convert rate to mbarn
tempNorm = 1/8.94e13
prob_observed(0.002117712, 3.286895e-05, 0.0005940328, 0.0001595911, 5.195151e-05,  0,  0, 0.0001977167, 0.005153151, 1.0, 1.0, 0.20)

# starting values n1, n2, eff....
fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik2/prob_dist_beam_100_Hendrik.RData"
load(file = fileName)
prob_dist_beam = prob_dist
# plot_ly(data = prob_dist_beam, x = ~n1, y = ~eff, color = ~prob, size = ~prob)
# plot_ly(data = prob_dist_beam, x = ~n1, y = ~n2, color = ~prob, size = ~prob) 

a1_mean = 40/factorCS
da1 =     20/factorCS
a2_mean = 780/factorCS
da2 =     50/factorCS

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
    a1 = runif(1, 0, a1_mean+da1)
    # a1 = 0
    prob_dist[1, "a1"] = a1
    a2 = runif(1, a2_mean-da2, a2_mean+da2)
    # a2 = 0
    prob_dist[1, "a2"] = a2
    effr = runif(1, 0.17, 0.27)
    #effr = 0.22
    prob_dist[1, "effr"] = effr
    
    prob_dist[1, "prob"] = prob_observed(n1, n2, b1, b2, b3, b4, b5, a1, a2, effn, effb, effr)
    # if (prob_dist[1, "prob"] > 1e-10) return ( prob_dist )
    return ( prob_dist )
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


a2a1 <- plot_ly(data = prob_dist, x = ~a1, y = ~a2, color = ~prob, size = ~prob, 
                text = ~paste0("Prob: ", format(prob,digits=2)),
                type = 'scatter', mode = 'markers')
a1effb <- plot_ly(data = prob_dist, x = ~effb, y = ~a1, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers')
a2effb <- plot_ly(data = prob_dist, x = ~effb, y = ~a2, color = ~prob, size = ~prob, 
                  text = ~paste0("Prob: ", format(prob,digits=2)),
                  type = 'scatter', mode = 'markers')
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
sb3 <- subplot( style(a2a1, showlegend = FALSE), style(a1effb, showlegend = FALSE), style(a2effb, showlegend = FALSE), 
        nrows = 3, margin = 0.01, titleX = T, titleY = T)
sb4 <- subplot( style(a2a1CS, showlegend = FALSE), style(a1effr, showlegend = FALSE), style(a2effr, showlegend = FALSE), 
                nrows = 3, margin = 0.01, titleX = T, titleY = T)

subplot(sb1, sb2,sb3, sb4, margin = 0.05, titleX = T, titleY = T)

# 3D plot ----------------
plot_ly(data = prob_dist, x = ~a1*factorCS, y = ~a2*factorCS, z = ~prob, color = ~prob, size = ~prob, 
        type="scatter3d", mode='markers' ) %>%
  layout(scene = list(xaxis = list(title = 'a1'),
                      yaxis = list(title = 'a2'),
                      zaxis = list(title = 'prob')))

fileName <- "~/Dropbox/Courses/R/alpha_n/Hendrik2/prob_dist_reaction_100_Hendrik.RData"
save(prob_dist, file = fileName)

# Testing ------------
maxP = prob_dist[prob_dist$prob == max(prob_dist$prob),]
prob_dist = prob_dist[prob_dist$prob > maxP$prob*0.0001,]




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

plot(0:1000, 0:1000, type = "n", xlim=c(0, 550), ylim=c(200, 950), xlab="a1n [mb]", ylab="a2n [mb]", axes = F)  # setting up coord. system
axis(side = 1, at = c(0,100,200,300,400,500))
axis(side = 2, at = c(200,300,400,500,600,700,800,900))
box()
points(x = maxP$a1*factorCS, y = maxP$a2*factorCS, col = "red")
points(x = 34, y = 780, col = "blue")
contour(x, y, cumDensity,
        xlab="a1n [mb]", ylab="a2n [mb]",
        levels=c(0.59, .80), lwd=4, lty=2, add = TRUE, col = "pink")


