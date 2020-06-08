### Sara M Schaal 
## Walsh & Lynch Eq. 5.6c

#### Model parameters ####
n <- c(2, 10)                          # number of loci 
enVar <- c(0.05, 0.1)                 # environmental variance added to the phenotype
ai <- c(0.01, 0.1)                     # effect size of the allele B 
pi <- seq(0, 1, by = 0.1)              # frequency of the allele B
m <- 1                                 # baseline value for mean
#########################

#### Build Dataframe ####
df.eq5.6c <- expand.grid("alpha" = ai, "freq" = pi, "enVar" = enVar, "m" = m, "nLoci" = n)
df.eq5.6c$traitMean <- NULL
df.eq5.6c$traitVar <- NULL

for(i in 1:nrow(df.eq5.6c)){
  df.eq5.6c$traitMean[i] <- df.eq5.6c$m[i] + 2*df.eq5.6c$alpha[i]*df.eq5.6c$freq[i]*df.eq5.6c$nLoci[i]
  df.eq5.6c$traitVar[i] <- 2*((df.eq5.6c$alpha[i])^2*df.eq5.6c$freq[i])*
    (1-df.eq5.6c$freq[i])*df.eq5.6c$nLoci[i] + df.eq5.6c$enVar[i]
}

#########################

#### Plotting ####
par(mfrow = c(2,2))
for(i in 1:length(ai)){
  ## first alpha value & number of loci##
  plot(traitMean~freq, data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                                          df.eq5.6c$enVar == enVar[1] & df.eq5.6c$nLoci == n[1],],
       col = "blue", type = "l", ylim = c(0.8,3.2), 
       main = paste("alpha = ", ai[i], "\nnLoci =", n[1]), 
       xlab = "Frequency of Allele", ylab = "Trait Mean")
  lines(traitMean+traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[1] & df.eq5.6c$nLoci == n[1],], 
        col = "blue", lty = 2)
  lines(traitMean-traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar ==enVar[1] & df.eq5.6c$nLoci == n[1],], 
        col = "blue", lty = 2)
  
  ## second environmental variance ##
  lines(traitMean~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[1],], 
        col = "green")
  lines(traitMean+traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[1],], 
        col = "green", lty = 2)
  lines(traitMean-traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[1],], 
        col = "green", lty = 2)
  if(i == 1){
    legend("topright", legend = c("env = 0.01",  "Var, env = 0.01", "env = 0.05", "Var, env = 0.05"), col = c("blue", "blue", "green", "green"), lty = c(1,2,1,2), cex = 0.7)
  }
  
  ## second alpha value & number of loci ###
  plot(traitMean~freq, 
       data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                          df.eq5.6c$enVar == enVar[1] & df.eq5.6c$nLoci == n[2],], 
       col = "blue", type = "l", ylim = c(0.8,3.2), 
       main = paste("alpha = ", ai[i], "\nnLoci =", n[2]), 
       xlab = "Frequency of Allele", ylab = "Trait Mean")
  
  lines(traitMean+traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[1] & df.eq5.6c$nLoci == n[2],], 
        col = "blue", lty = 2)
  lines(traitMean-traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[1] & df.eq5.6c$nLoci == n[2],], 
        col = "blue", lty = 2)
  
  lines(traitMean~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[2],], 
        col = "green")
  lines(traitMean+traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[2],], 
        col = "green", lty = 2)
  lines(traitMean-traitVar~freq, 
        data = df.eq5.6c[df.eq5.6c$alpha == ai[i] & 
                           df.eq5.6c$enVar == enVar[2] & df.eq5.6c$nLoci == n[2],], 
        col = "green", lty = 2)
  
  
}
