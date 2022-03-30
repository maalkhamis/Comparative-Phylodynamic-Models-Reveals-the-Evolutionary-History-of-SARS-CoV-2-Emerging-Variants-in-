

#------------------------------------------------------------------------
#Code to run Skygrowth for SARS-coV2 genomes
#------------------------------------------------------------------------

#Written by Nick Fountain-Jones (Nick.FountainJones@utas.edu.au) with help from Xavier Didelot and Erik VOlz

library(treeio)
library(ggtree)
library(tidyverse)
library(treestructure) 
library(phylodyn)# needs to options(buildtools.check = function(action) TRUE ) before installing from Git
#install_github("mdkarcher/phylodyn")
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(skygrowth)
library(coda)
library(ape)
library(treedater)
library(lubridate)
library(ggplot2)

####################################################################
#Kappa
####################################################################

BEASTtreKappa <-read.nexus('Kappa_Ex_UCED_Sym2.tre')

ggtree(BEASTtreKappa ) + geom_tiplab(size=3)
#check for any clades of interest within each MCC tree

#mcmc provides more accurate CIs
Kappagrowth <- skygrowth.mcmc(BEASTtreKappa,  res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 3e+07, control=list(thin=1e3) ) 

#check convergence
KapagrowthMCMC <- as.mcmc(cbind(Kappagrowth $growthrate[,1:(ncol(Kappagrowth $growthrate)-1)],Kappagrowth$ne,Kappagrowth $tau))
effectiveSize(KapagrowthMCMC)

#
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(Kappagrowth  )+theme_bw()
neplot(Kappagrowth)+theme_bw()
R.plot(Kappagrowth , forward=TRUE, gamma=0.90)+theme_bw()

#save files
save(Kappagrowth, file="Kappa")
load("global_covid")

#compare to phylodyn
GlobalBSP_kappa <- BNPR(BEASTtreKappa, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP )

####################################################################
#Delta
####################################################################

#need to update these - negative branch lengths on the MCC

BEASTtreDelta <-read.nexus('Delta_Eg_UCED_ASym1_caHeights')

ggtree(BEASTtreDelta  ) + geom_tiplab(size=3)
#check for any clades of interest within each MCC tree

#mcmc provides more accurate CIs
Deltagrowth <- skygrowth.mcmc(BEASTtreDelta ,  res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
DeltagrowthMCMC <- as.mcmc(cbind(Deltagrowth $growthrate[,1:(ncol(Deltagrowth$growthrate)-1)],Deltagrowth$ne,Deltagrowth$tau))
effectiveSize(DeltagrowthMCMC)

#
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(Deltagrowth)+theme_bw()
neplot(Deltagrowth)+theme_bw()
R.plot(Deltagrowth, forward=TRUE, gamma=0.90)+theme_bw()

#save files
save(Deltagrowth, file="Delta")
load("Delta")

#compare to phylodyn
GlobalBSP_delta <- BNPR(BEASTtreDelta, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP_delta )


####################################################################
#Alpha
####################################################################
BEASTtreAlpha <-read.nexus('Alpha_Eg_UCED_ASym2_caHeights')

ggtree(BEASTtreAlpha) + geom_tiplab(size=3)
#check for any clades of interest within each MCC tree

#mcmc provides more accurate CIs
Alphagrowth <- skygrowth.mcmc(BEASTtreAlpha ,  res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
AlphagrowthMCMC <- as.mcmc(cbind(Alphagrowth$growthrate[,1:(ncol(Alphagrowth$growthrate)-1)],Alphagrowth$ne,Alphagrowth$tau))
effectiveSize(AlphagrowthMCMC)

#
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(Alphagrowth)+theme_bw()
neplot(Alphagrowth)+theme_bw()
R.plot(Alphagrowth, forward=TRUE, gamma=0.90)+theme_bw()

#compare to phylodyn
GlobalBSP_alpha <- BNPR(BEASTtreAlpha, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                        beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                        derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP_alpha )

save(Alphagrowth, file="Alpha")
load("global_covid")

####################################################################
#Beta
####################################################################
BEASTtreBeta <-read.nexus('Beta_Eg_UCED_Sym2.tre')

ggtree(BEASTtreBeta) + geom_tiplab(size=3)
#check for any clades of interest within each MCC tree

#mcmc provides more accurate CIs
Betagrowth <- skygrowth.mcmc(BEASTtreAlpha ,  res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
BetagrowthMCMC <- as.mcmc(cbind(Alphagrowth$growthrate[,1:(ncol(Alphagrowth$growthrate)-1)],Alphagrowth$ne,Alphagrowth$tau))
effectiveSize(AlphagrowthMCMC)

#plots

growth.plot(Betagrowth)+theme_bw()
neplot(Betagrowth)+theme_bw()
R.plot(Betagrowth, forward=TRUE, gamma=0.90)+theme_bw()

#compare to phylodyn
GlobalBSP_beta <- BNPR(BEASTtreBeta, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                        beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                        derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP_beta )

save(Betagrowth, file="Beta")

####################################################################
#Eta
####################################################################


BEASTtreEta <- read.nexus('Eta_Eg_UCED_Sym1+updated.tre')


ggtree(BEASTtreEta ) + geom_tiplab(size=3)

EtaGrowth <-  skygrowth.mcmc(BEASTtreEta,  res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 3e+07, control=list(thin=1e3) ) 


#check convergence
EtagrowthMCMC <- as.mcmc(cbind(EtaGrowth $growthrate[,1:(ncol(EtaGrowth $growthrate)-1)],EtaGrowth$ne,Kappagrowth $tau))
effectiveSize(EtagrowthMCMC )

#
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(EtaGrowth   )+theme_bw()
neplot(EtaGrowth)+theme_bw()
R.plot(EtaGrowth , forward=TRUE, gamma=0.90)+theme_bw()

#save files
save(TwentyDgrowth, file="20D")
load("global_covid")

#compare to phylodyn
GlobalBSP_eta <- BNPR(BEASTtreEta, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP_eta )

#Previous analyses
####################################################################
#Lineage 20G
####################################################################
BEASTtre20G <-read.nexus('20GAliSym2.tre')
ggtree(BEASTtre20G ) + geom_tiplab(size=3)
#check for any clades of interest within each MCC tree
#treeStBEAST <-  trestruct(BEASTtre , minCladeSize = 100, minOverlap = -Inf, nsim = 10000,
                         # level = 0.05, ncpu = 1, verbosity = 1) 


#plot(treeStBEAST, use_ggtree = TRUE) 

TwentyGsg<- skygrowth.map(BEASTtre20G, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T)) #res=35

plot(TwentyDsg)

#mcmc provides more accurate CIs
TwentyGgrowth <- skygrowth.mcmc(BEASTtre20G, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
TwentyGgrowthMCMC <- as.mcmc(cbind(TwentyDgrowth$growthrate[,1:(ncol(TwentyDgrowth$growthrate)-1)],TwentyDgrowth$ne,TwentyDgrowth$tau))
effectiveSize(TwentyGgrowthMCMC)

# need date info
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(TwentyGgrowth)+theme_bw()
neplot(TwentyGgrowth)+theme_bw()
R.plot(TwentyGgrowth, forward=TRUE, gamma=0.90)+theme_bw()

#save files
save(TwentyGgrowth, file="20G")
load("global_covid")

#compare to phylodyn
GlobalBSP <- BNPR(BEASTtre20G, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP )

####################################################################
#Lineage V1
####################################################################
BEASTtreV1 <-read.nexus('V1aliASym1.tre')
ggtree(BEASTtreV1 ) + geom_tiplab(size=3) #neg branch lengths

#check for any clades of interest within each MCC tree
#treeStBEAST <-  trestruct(BEASTtre , minCladeSize = 100, minOverlap = -Inf, nsim = 10000,
# level = 0.05, ncpu = 1, verbosity = 1) 


#plot(treeStBEAST, use_ggtree = TRUE) 

V1sg<- skygrowth.map(BEASTtreV1, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T)) #res=35

plot(TwentyDsg)

#mcmc provides more accurate CIs
V1growth <- skygrowth.mcmc(BEASTtreV1, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
V1growthMCMC <- as.mcmc(cbind(V1growth$growthrate[,1:(ncol(V1growth$growthrate)-1)],V1growth$ne,V1growth $tau))
effectiveSize(V1growthMCMC )

# need date info
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(V1growth)+theme_bw()
neplot(V1growth)+theme_bw()
R.plot(V1growth, forward=TRUE, gamma=0.90)+theme_bw()

#save files
save(TwentyGgrowth, file="20G")
load("global_covid")

#compare to phylodyn
GlobalBSP <- BNPR(BEASTtreV1, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP )

####################################################################
#Lineage V2
####################################################################
BEASTtreV2 <-read.nexus('V2aliSym1.tre')
ggtree(BEASTtreV2 ) + geom_tiplab(size=3) #neg branch lengths

#check for any clades of interest within each MCC tree
#treeStBEAST <-  trestruct(BEASTtre , minCladeSize = 100, minOverlap = -Inf, nsim = 10000,
# level = 0.05, ncpu = 1, verbosity = 1) 


#plot(treeStBEAST, use_ggtree = TRUE) 

V2sg<- skygrowth.map(BEASTtreV2, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T)) #res=35

plot(TwentyDsg)

#mcmc provides more accurate CIs
V2growth <- skygrowth.mcmc(BEASTtreV2, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

#check convergence
V2growthMCMC <- as.mcmc(cbind(V2growth$growthrate[,1:(ncol(V1growth$growthrate)-1)],V1growth$ne,V1growth $tau))
effectiveSize(V2growthMCMC )

# need date info
sampleDates <- treedater::sampleYearsFromLabels( ebov_ml$tip.label, delimiter='|' )

#plots

growth.plot(V2growth)+theme_bw()
neplot(V2growth)+theme_bw()
R.plot(V2growth, forward=TRUE, gamma=0.90)+theme_bw()

computeR(V2growth, gamma=0.90)

#save files
save(V2growth, file="V2")
load("global_covid")

#compare to phylodyn
GlobalBSP <- BNPR(BEASTtreV2, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                  beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                  derivative = FALSE, forward = TRUE)
plot_BNPR(GlobalBSP )

