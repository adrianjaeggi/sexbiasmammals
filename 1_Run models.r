### load packages
library(rlang)
library(Rcpp)
library(rstan)
library(rethinking)
library(brms)

#stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

memory.limit(size=100000)

## read data and phylo sample
d<- read.csv("Sex bias mammals data.csv")

trees<- read.tree("Sex bias mammals 100 phylos.tre")

## note the following changes in species names to match our original database to the vertlife phylo
# Capra_aegagrus_hircus --> Capra_aegagrus
# Stenella_longirostris_longirostris --> Stenella_longirostris
# Tursiops_spp --> Tursiops_truncatus
# Panthera_leo_leo & Panthera_leo_persica --> Panthera_leo
# Helogale_paruvla --> Helogale_parvula (fix in google doc since clearly a typo)
# Equus_ferus_przewalskii --> Equus_ferus
# Equus_burchellii_quagga --> Equus_quagga
# Cercopithecus_kandti & Cercopithecus_mitis_erythrarchus & Cercopithecus_mitis_stuhlmanni --> Cercopithecus_mitis
# Piliocolobus_kirkii --> Procolobus_kirkii
# Gorilla_beringei_beringei --> Gorilla_beringei
# Gorilla_gorilla_gorilla --> Gorilla_gorilla
# Eulemur_fulvus_rufus --> Eulemur_fulvus
# Hapalemur_griseus_alaotensis --> Hapalemur_griseus

# convert to covariance matrix (see https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html)
A<- list() 
for(i in 1:100){
  A[[i]]<- vcv.phylo(trees[[i]])
}

# duplicate species list for mixed-effects model
d$phylo<- d$Genus_species

## set prior
prior.sexbias<- c(
  prior(normal(0,5), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,5), "Intercept", dpar="muMalebiased"),
  prior(normal(0,2), "b", dpar="muFemalebiased"), 
  prior(normal(0,2), "b", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

### Model 1: between-group conflict
bgc.1.loop1<- brm(Between.Conflict_zero_center.fact ~ 1 + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# check summary, plot for diagnostics --> all good

# extract samples
post_bgc.1_looped<- posterior_samples(bgc.1.loop1)

# loop over the rest of the tree sample
for(i in 2:100){
  bgc.1.loopi<- brm(Between.Conflict_zero_center.fact ~ 1 + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
  post_bgc.1.loopi<- posterior_samples(bgc.1.loopi)
  post_bgc.1_looped<- rbind(post_bgc.1_looped, post_bgc.1.loopi)
  save(post_bgc.1_looped, file="post_bgc.1_looped.robj")
}

nrow(post_bgc.1_looped)/4000 # check all 100 models added 


### Model 1: Movement
mov.1.loop1<- brm(Movement_zero_center.fact ~ 1 + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
post_mov.1_looped<- posterior_samples(mov.1.loop1)
for(i in 2:100){
  mov.1.loopi<- brm(Movement_zero_center.fact ~ 1 + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
  post_mov.1.loopi<- posterior_samples(mov.1.loopi)
  post_mov.1_looped<- rbind(post_mov.1_looped, post_mov.1.loopi)
  save(post_mov.1_looped, file="post_mov.1_looped.robj")
}
nrow(post_mov.1_looped)/4000 # check all 100 models added 

## prep for Model 2: center and re-scale dimorphism
d$Sex_Dim.z<- (d$Sex_Dim-1) / sd(d$Sex_Dim) # 0 = monomorphic, units = SD


### Model 2: Between-group conflict
bgc.2.loop1<- brm(Between.Conflict_zero_center.fact ~ Sex_Dim.z + Food_resource_defendable + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 5000, warmup = 4000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.0001), inits=0)
post_bgc.2_looped<- posterior_samples(bgc.2.loop1)
for(i in 2:100){
  bgc.2.loopi<- brm(Between.Conflict_zero_center.fact ~ Sex_Dim.z + Food_resource_defendable + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 5000, warmup = 4000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.0001), inits=0)
  post_bgc.2.loopi<- posterior_samples(bgc.2.loopi)
  post_bgc.2_looped<- rbind(post_bgc.2_looped, post_bgc.2.loopi)
  save(post_bgc.2_looped, file="post_bgc.2_looped.robj")
}
nrow(post_bgc.2_looped)/4000 # check all 100 models added

### Model 2: Movement
mov.2.loop1<- brm(Movement_zero_center.fact ~ Sex_Dim.z + Food_resource_defendable + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20, stepsize = 0.00001), inits=0)
post_mov.2_looped<- posterior_samples(mov.2.loop1)
for(i in 2:100){
  mov.2.loopi<- brm(Movement_zero_center.fact ~ Sex_Dim.z + Food_resource_defendable + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20, stepsize = 0.00001), inits=0)
  post_mov.2.loopi<- posterior_samples(mov.2.loopi)
  post_mov.2_looped<- rbind(post_mov.2_looped, post_mov.2.loopi)
  save(post_mov.2_looped, file="post_mov.2_looped.robj")
}
nrow(post_mov.2_looped)/4000 # check all 100 models added 


### Model 3: Between-group conflict
bgc.3.loop1<- brm(Between.Conflict_zero_center.fact ~ Primate + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.00001), inits=0)
post_bgc.3_looped<- posterior_samples(bgc.3.loop1)
for(i in 2:100){
  bgc.3.loopi<- brm(Between.Conflict_zero_center.fact ~ Primate + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.00001), inits=0)
  post_bgc.3.loopi<- posterior_samples(bgc.3.loopi)
  post_bgc.3_looped<- rbind(post_bgc.3_looped, post_bgc.3.loopi)
  save(post_bgc.3_looped, file="post_bgc.3_looped.robj")
}
nrow(post_bgc.3_looped)/4000 # check all 100 models added 

### Model 3: Movement
mov.3.loop1<- brm(Movement_zero_center.fact ~ Primate + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20, stepsize = 0.0001), inits=0)
post_mov.3_looped<- posterior_samples(mov.3.loop1)
for(i in 2:100){
  mov.3.loopi<- brm(Movement_zero_center.fact ~ Primate + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20, stepsize = 0.0001), inits=0)
  post_mov.3.loopi<- posterior_samples(mov.3.loopi)
  post_mov.3_looped<- rbind(post_mov.3_looped, post_mov.3.loopi)
  save(post_mov.3_looped, file="post_mov.3_looped.robj")
}
nrow(post_mov.3_looped)/4000 # check all 100 models added 


### Model 4
bgc.4.loop1<- brm(Between.Conflict_zero_center.fact ~ Movement_zero_center.fact + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[1]]), 
                  prior = prior.sexbias, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.0001), inits=0)
post_bgc.4_looped<- posterior_samples(bgc.4.loop1)
for(i in 2:100){
  bgc.4.loopi<- brm(Between.Conflict_zero_center.fact ~ Movement_zero_center.fact + (1|Genus_species) + (1|phylo), data = d, family = categorical(), cov_ranef = list(phylo = A[[i]]), 
                    prior = prior.sexbias, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 19, stepsize = 0.0001), inits=0)
  post_bgc.4.loopi<- posterior_samples(bgc.4.loopi)
  post_bgc.4_looped<- rbind(post_bgc.4_looped, post_bgc.4.loopi)
  save(post_bgc.4_looped, file="post_bgc.4_looped.robj")
}
nrow(post_bgc.4_looped)/12000 # check all 100 models added (note larger n samples per model here)


