 
#### Excercise 5.7 (page 12X) ####
# Calculating the marginal effects (W_AA) of 
# homozygous A (AA) individuals, while considering the 
# impact of gametic-phase disequilibrium (i.e. linkage disequilibrium)

## Before we get started a quick refresher
# on conditional probabilities.

# Consider Pr(AABB|AA) - 
# This is equivalent to the probability of the two-locus
# genotype AABB given the probability of AA.
# This can also be obtained with the formula Pr(AABB,BB)/Pr(BB)
# Which reads as the probability of AABB AND BB divided by the 
# probability of BB.

## Libraries 
library(ggplot2)
library(cowplot)

## Starting parameters
# Allele frequencies of dominant allele at each locus (A or B)
p = 0.5   # Freq. A
q = 0.5   # Freq. B
## Probability of A and B together
pq = 1 #p*q

# Linkage Disequilibrium between loci
(LD = pq - (p*q))
# Calculated as the probability of haplotype AB - probabilty of allele A x pr. of allele B

# Fitness of various genotypes (assumes A locus is AA)
W_BB <- 1
W_Bb <- 1
W_bb <- 1

marg_fitness <- function(p=0.5,q=0.5,LD=0.5,
                         W1=1,W2=1,W3=1){
y <-  W1*((p*q + LD)^2/p^2) + 
      W2*((2*(p*q + LD)*(p*(1-q)- LD))/p^2) +
      W3*((p*(1-q) - LD)^2/p^2)
return(y)
}

marg_fitness(p,q,LD,W_BB,W_Bb,W_bb)

#### Exploring the effect of LD on marginal fitness
# Equal loci 2 fitnesses (W_BB = W_Bb = W_bb = 1)
LD = seq(0,1,by=0.1)
plot(marg_fitness(p,q,LD,W_BB,W_Bb,W_bb)~LD,xlab="Linkage Disequilibrium",ylab="W_AA")
# Linkage DisEq. doesn't impact marginal fitness of AA IF fitness of all genotype combinations at loci 2 are the same

## Unequal loci fitnesses (W_BB = W_B_b = 2, W_bb = 1)
plot(marg_fitness(p,q,LD,2,2,1)~LD,xlab="Linkage Disequilibrium",ylab="W_AA")
## Unequal loci fitnesses (W_BB = 1, W_B_b = 2, W_bb = 1)
plot(marg_fitness(p,q,LD,1,2,1)~LD,xlab="Linkage Disequilibrium",ylab="W_AA")
## Unequal loci fitnesses (W_BB = 1, W_B_b = 1, W_bb = 2)
plot(marg_fitness(p,q,LD,1,1,2)~LD,xlab="Linkage Disequilibrium",ylab="W_AA")


### Marginal fitness based on variable starting allele frequencies (p) and
  # several different levels of LD
LD = c(0.1,0.5,1)
p = seq(0.1,1,by=0.01)
q = 0.8

out <- NULL
for(i in 1:length(LD)) out <- rbind(out,cbind(p=p,LD=LD[i],fitness=marg_fitness(p,q,LD[i],2,2,1)))
out <- data.frame(out) 

ggplot(out,aes(x=p,y=fitness,colour=as.factor(LD))) + 
  geom_point() +
  theme_cowplot() +
  labs(colour="LD")

## LD really only impact marginal fitness when starting allele frequencies
 # are small
