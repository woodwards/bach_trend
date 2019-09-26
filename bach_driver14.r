# translate runDREAMZS_hmixbatch.m into R
# Simon Woodward, DairyNZ 2017-2018
# now allows keeping previous results or starting from scratch
# save runlist including diagnostics
# bach_driver11.r + bach_script11.stan: 
# - apply priors to parameters (requires intial values to chains)
# - unique name for log file
# - increase weighting of llflow
# - reorganise parameter trace plots

# remove all variables
rm(list=ls()) 
# dev.off() # close figures if they exist

suppressMessages({
  library(rstan)
  library(tidyverse)
  library(cowplot)
  library(parallel)
  library(installr)
  library(truncnorm)
  library(ggthemes)
})

# functions
"%notin%" <- function(x,y)!("%in%"(x,y))

# which run directory?
out_path <- 'run - new_trend/'
print(out_path)

# kill any processes that didn't terminate
kill_all_Rscript_s()

# open log_file
log_file_name <- paste(out_path, 'log', format(Sys.time(), "_%Y%m%d_%H%M"), '.txt', sep='')
log_file <- file(log_file_name, open='w') # delete duplicate log file if necessary
close(log_file)
log_file <- file(log_file_name, open='a') # open file for logging

# model and parameters
stan_model <- paste(out_path, 'bach_script13.stan', sep='')
stan_pars_raw <- c('medb0raw', 'medd1raw', 'slowb0raw', 'slowd1raw',
                   'chem1fastraw0', 'chem1medraw0', 'chem1slowraw0', 
                   'chem2fastraw0', 'chem2medraw0', 'chem2slowraw0',
                   'chem1fastraw1', 'chem1medraw1', 'chem1slowraw1', 
                   'chem2fastraw1', 'chem2medraw1', 'chem2slowraw1')
stan_pars <- c('medb0', 'meda1', 'slowb0', 'slowa1', 
               'medBFImax', 'medrec', 'slowBFImax', 'slowrec', 
               'chem1fast0', 'chem1med0', 'chem1slow0', 
               'chem2fast0', 'chem2med0', 'chem2slow0',
               'chem1fast1', 'chem1med1', 'chem1slow1', 
               'chem2fast1', 'chem2med1', 'chem2slow1',
               'llchem1', 'llchem2')

# record version information
cat(R.version.string, file=log_file, sep='\n')
cat(paste('Stan version', stan_version()), file=log_file, sep='\n')

# init rstan
rstan_options(auto_write = TRUE)
options(mc.cores = detectCores())

# read data files
source('read_data9.r')
cat(paste(nrow(runlist),'lines in runlist.dat'), file=log_file, sep='\n')

# prepare output files
nruns <- nrow(runlist)
from_scratch <- all(runlist$control %notin% "keep")
source('prep_output8.r') 

# loop through runlist
rows <- 1:nruns # to run all sites
# rows <- c(12) # to run a subset of sites, useful for testing
i <- rows[[1]] # set i useful for line by line testing
for (i in rows) {
  
  # seed
  set.seed(i)
  
  # write header
  data_file_name <- paste(runlist$catchfile[i], '_data.dat', sep='')
  opt_file_name <- runlist$optfile[i]
  header <- paste(i, data_file_name, opt_file_name, format(Sys.time(), "%a %d %b %H:%M:%S %Y"))
  print(header)
  cat(header, file=log_file, sep='\n')
  
  # get mode
  mode <- runlist[i, 'control']
  if (mode=='zero') {
    
    # report
    cat(paste('max(Rhat)  = zero'), file=log_file, sep='\n')
    cat(paste('elapsed(h) = 0'), file=log_file, sep='\n')
    
    # insert into output files
    irows <- (i*7-6):(i*7)
    bestxs[i, ] <- 0
    quartxs[irows, ] <- 0
    sbestxs[i, ] <- 0
    squartxs[irows, ] <- 0
    
  } else if (mode=='calib') {
    
    # assemble run options/data into vectors (?)
    arun <- runlist[i, ]
    adata <- data[grep(pattern=data_file_name, x=data$file), ]
    aarea <- tibble(area=adata$area[1])
    aoptions <- options[opt_file_name, ]
    aalloptions <- cbind(arun, aoptions, aarea) # combine into one data table
    startcalib <- aalloptions$startcalib
    endcalib <- aalloptions$endcalib
    catchname <- aalloptions$catchname
    setname <- aalloptions$setname
    keeps <- c('startrun', 'startcalib', 'endcalib', 'startvalid', 'endvalid')
    intoptions <- as.vector(t(aalloptions[keeps]))
    nintoptions <- length(intoptions)
    keeps <- c('chem1ae', 'chem1re', 'chem2ae', 'chem2re', 'area')
    realoptions <- as.vector(t(aalloptions[keeps]))
    nrealoptions <- length(realoptions)
    date <- as.vector(adata$date)
    ndate <- length(date)
    flow <- as.vector(adata$flow)
    TP <- as.vector(adata$TP)
    TN <- as.vector(adata$TN)
    # rain <- as.vector(adata$rain)
    # pet <- as.vector(adata$pet)
    meanTPcalib <- mean(TP[startcalib:endcalib], na.rm=TRUE)
    meanTNcalib <- mean(TN[startcalib:endcalib], na.rm=TRUE)
    TP[is.na(TP)] <- -1 # stan doesn't understand NA
    TN[is.na(TN)] <- -1 # stan doesn't understand NA
    
    # # catch bad pet
    # j <- which(pet < 0)
    # pet[j] <- (pet[j-1] + pet[j+1]) / 2
    # stopifnot(all(pet>=0))
    
    # 
    nits <- aalloptions[1, 'nits']
    
    # fit the model
    wits <- nits # warmup iterations (minimum recommended = 150)
    sits <- 200 # sampling iterations FIXME enough? was 400 for bach
    
    # // priors
    priorwide <- 0.3
    priornarrow <- 0.1
    # medb0raw ~ normal(1.0, priornarrow); 
    # medd1raw ~ normal(0.4, priorwide); 
    # slowb0raw ~ normal(1.0, priornarrow); 
    # slowd1raw ~ normal(0.7, priorwide); 
    # chem1fastraw ~ normal(0.2/2.0, priorwide); 
    # chem1medraw ~ normal(0.1/2.0, priorwide); 
    # chem1slowraw ~ normal(0.1/2.0, priorwide); 
    # chem2fastraw ~ normal(2.0/12.0, priorwide);
    # chem2medraw ~ normal(3.0/12.0, priorwide); 
    # chem2slowraw ~ normal(1.0/12.0, priorwide);
    
    chain_init <- function(){list(
      medb0raw=rtruncnorm(1, 0, 1, 1.0, priornarrow),
      medd1raw=rtruncnorm(1, 0, 1, 0.4, priorwide),
      slowb0raw=rtruncnorm(1, 0, 1, 1.0, priornarrow),
      slowd1raw=rtruncnorm(1, 0, 1, 0.7, priorwide),
      chem1fastraw0=rtruncnorm(1, 0, 1, 0.2/2.0, priorwide),
      chem1medraw0=rtruncnorm(1, 0, 1, 0.1/2.0, priorwide),
      chem1slowraw0=rtruncnorm(1, 0, 1, 0.1/2.0, priorwide),
      chem2fastraw0=rtruncnorm(1, 0, 1, 2.0/12.0, priorwide),
      chem2medraw0=rtruncnorm(1, 0, 1, 3.0/12.0, priorwide),
      chem2slowraw0=rtruncnorm(1, 0, 1, 1.0/12.0, priorwide),
      chem1fastraw1=rtruncnorm(1, 0, 1, 0.2/2.0, priorwide),
      chem1medraw1=rtruncnorm(1, 0, 1, 0.1/2.0, priorwide),
      chem1slowraw1=rtruncnorm(1, 0, 1, 0.1/2.0, priorwide),
      chem2fastraw1=rtruncnorm(1, 0, 1, 2.0/12.0, priorwide),
      chem2medraw1=rtruncnorm(1, 0, 1, 3.0/12.0, priorwide),
      chem2slowraw1=rtruncnorm(1, 0, 1, 1.0/12.0, priorwide)
    )}
    
    stan_init <- list(chain_init(), chain_init(), chain_init(), chain_init())
    
    fit <- stan(stan_model,
                data=c('nintoptions', 'intoptions', 'nrealoptions', 'realoptions',
                       'ndate', 'date', 'flow', 'TP', 'TN', 'meanTPcalib', 'meanTNcalib'),
                # control=list(adapt_delta=0.90, stepsize=0.005, max_treedepth=15),
                warmup=wits, 
                iter=wits+sits, 
                chains=4, 
                init=stan_init, # needed if priors
                # cores = 1, # disable parallel processing
                # save_warmup=FALSE, # warmup iterations required for some diagnostics
                seed=i)

    # memory management if necessary
    # save.image(file="temp.RData") 
    # rm(list=ls())
    # load(file="temp.RData")
    
    # extract samples
    samples <- as.data.frame(rstan::extract(fit)) %>% # very large!
      select(-starts_with('bach')) %>% # discard model traces. This is a bit slow.
      mutate(setname=setname) # add identifier 
    write_rds(samples, paste(out_path, setname, '_samples.rds', sep=''))    
    
    # extract quartiles
    results <- summary(fit)$summary # quartiles
    keeps <- grep('bach*', rownames(results), invert=FALSE) # extract model traces 'bach[]'
    traces <- as.data.frame(results[keeps,]) # don't use tibble because want row names
    keeps <- grep('bach*', rownames(results), invert=TRUE) # discard model traces 'bach[]'
    results <- as.data.frame(results[keeps,]) # don't use tibble because want row names
    maxRhat <- max(results$Rhat, na.rm=TRUE) # should be < 1.2 
    warning <- ''
    if (maxRhat >= 1.1) { # generate a warning and save fit
      warning <- paste('  Warning: max(Rhat) =', maxRhat, '>= 1.l')
      print(warning)
      # http://discourse.mc-stan.org/t/saving-and-sharing-an-rstan-model-fit/1059
      # fit@stanmodel@dso <- new("cxxdso") # have to kill the Dynamic Shared Object inside the stanmodel slot? 
      # write_rds(fit, paste(out_path, setname, '_fit.rds', sep='')) # big file!
    }
    runtime <- as_tibble(get_elapsed_time(fit)) %>% mutate(total=warmup+sample)
    resultst <- cbind.data.frame(nms=names(results), t(results), stringsAsFactors=FALSE) # transpose, keeping names
    # don't use mutate; it kills row names
    resultst[, 1] <- 0 # replace first columns 'nms' with 0
    resultst[, 'totalgens'] <- wits + sits
    resultst[, 'gelmanr'] <- maxRhat
    resultst[, 'elapsed'] <- max(runtime$total)/3600 # hours
    resultst[, 'setname'] <- setname
    
    write_rds(traces, paste(out_path, setname, '_traces.rds', sep=''))
    write_rds(resultst, paste(out_path, setname, '_quartiles.rds', sep=''))
    
    
    # report
    cat(paste('iterations =', wits+sits), file=log_file, sep='\n')
    cat(paste('max(Rhat)  =', maxRhat, warning), file=log_file, sep='\n')
    cat(paste('elapsed(h) =', max(runtime$total)/3600), file=log_file, sep='\n')
    runrecord$gelmanr[i] <- maxRhat
    runrecord$elapsed[i] <- max(runtime$total)/3600
    
    # insert into output files
    irows <- (i*7-6):(i*7)
    keeps <- match(parcols, colnames(resultst), nomatch=1) # col numbers we want, use col 1 if missing
    bestxs[i, ] <- resultst['50%', keeps]
    quartxs[irows, ] <- resultst[c('2.5%','2.5%','25%','50%','75%','97.5%','97.5%'), keeps]
    keeps <- match(statcols, colnames(resultst), nomatch=1) # col numbers we want, use col 1 if missing
    sbestxs[i, ] <- resultst['50%', keeps]
    squartxs[irows, ] <- resultst[c('2.5%','2.5%','25%','50%','75%','97.5%','97.5%'), keeps]
    
  } else { # do nothing
    
    # report
    cat(paste('iterations = 0'), file=log_file, sep='\n')
    cat(paste('max(Rhat)  = unchanged'), file=log_file, sep='\n')
    cat(paste('elapsed(h) = 0'), file=log_file, sep='\n')
    
  } # fit the model
  
  # write results so far
  write_tsv(bestxs, file.path(out_path, 'bestxs.tsv'), col_names=FALSE, na='#N/A') 
  write_tsv(quartxs, file.path(out_path, 'quartxs.tsv'), col_names=FALSE, na='#N/A') 
  write_tsv(sbestxs, file.path(out_path, 'sbestxs.tsv'), col_names=FALSE, na='#N/A') 
  write_tsv(squartxs, file.path(out_path, 'squartxs.tsv'), col_names=FALSE, na='#N/A') 
  write_tsv(runrecord, file.path(out_path, 'runrecord.tsv'), col_names=TRUE) 
  
  # continue
  i <- i + 1  
  
} # end loop

#
cat(format(Sys.time(), "%a %d %b %H:%M:%S %Y"), file=log_file, sep='\n')
cat("\n")
print("Stan compile information:")
print(fit@stanmodel@dso) # print compilation information
close(log_file)

#
showConnections(all=FALSE)
closeAllConnections() 
print("Finished. Generating diagnostics and charts for last run.")

# save last image
# save.image(file="temp.RData")
# rm(list=ls())
# load(file="temp.RData")

# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html 
cat("\n")
# get_adaptation_info(fit)
print(fit, stan_pars_raw)
stan_diag(fit)
stan_par(fit, stan_pars_raw[[2]]) # look at misbehaving parameters
cat("\n")

# https://mvuorre.github.io/post/2015/08-21-ggplot2-visualizing-prior-posterior/
# stan_hist(fit, stan_pars_raw)
stan_pars_calc <- c("medb0", "meda1", "slowb0", "slowa1",
                    "medBFImax", "medrec", "slowBFImax", "slowrec")
stan_pars_raw_plus <- c(stan_pars_raw, stan_pars_calc)
posteriors <- select(samples, stan_pars_raw_plus) %>%
  gather() %>% 
  mutate(dist="post",
         key=factor(key, levels=stan_pars_raw_plus))
npriorsamples <- 800
priors <- vector("list", npriorsamples)
for (i in 1:npriorsamples){
  priors[[i]] <- as_tibble(chain_init())
}
priors <- bind_rows(priors)
priors <- priors %>%
  mutate(medd1scale = (-0.1*log(1-0.5)) + ((-0.1*log(1-0.99)) - (-0.1*log(1-0.5))) * medd1raw,
         meda1 = (1-exp(-abs(medd1scale)*10)),
         medratio = 1/(1-meda1),
         medb0 = if_else(medratio>1 , medb0raw/medratio , medb0raw),
         medBFImax = medb0/(1-meda1),
         medrec = if_else(medb0<1 , meda1/(1-medb0) , 0),
         slowd1scale = (-0.1*log(1-0.99)) + ((-0.1*log(1-0.9999)) - (-0.1*log(1-0.99))) * slowd1raw,
         slowa1 = (1-exp(-abs(slowd1scale)*10)),
         slowratio = 1/(1-slowa1),
         slowb0 = if_else(slowratio>1 , slowb0raw/slowratio , slowb0raw),
         slowBFImax = slowb0/(1-slowa1),
         slowrec = if_else(slowb0<1 , slowa1/(1-slowb0) , 0)
  )
priors <- priors %>%
  select(stan_pars_raw_plus) %>%
  gather(priors) %>% 
  mutate(dist="prior") %>% 
  rename(key=priors) %>%
  mutate(key=factor(key, levels=stan_pars_raw_plus))
# dist <- rbind(posteriors, priors) 
# https://stackoverflow.com/questions/6957549/overlaying-histograms-with-ggplot2-in-r
file_name <- paste(out_path, 'last_post_prior_raw.png', sep="")
png(file_name, height=9, width=9, units='in', res=300)
p1 <- ggplot() +
  theme_solarized() +
  labs(title=setname) +
  geom_histogram(data=filter(priors, key %in% stan_pars_raw), mapping=aes(value), fill="red", alpha=0.5, boundary=1, binwidth=0.05) +
  geom_histogram(data=filter(posteriors, key %in% stan_pars_raw), mapping=aes(value), fill="blue", alpha=0.5, boundary=1, binwidth=0.05) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  facet_wrap("key", nrow=2)
print(p1)
dev.off()
file_name <- paste(out_path, 'last_post_prior_calc.png', sep="")
png(file_name, height=9, width=9, units='in', res=300)
p2 <- ggplot() +
  theme_solarized() +
  labs(title=setname) +
  geom_histogram(data=filter(priors, key %in% stan_pars_calc), mapping=aes(value), fill="red", alpha=0.5, boundary=1) +
  geom_histogram(data=filter(posteriors, key %in% stan_pars_calc), mapping=aes(value), fill="blue", alpha=0.5, boundary=1) +
  # scale_x_continuous(breaks=c(0,0.5,1)) +
  facet_wrap("key", nrow=2, scales="free")
print(p2)
dev.off()

# use modified code
source('pairs.stanfit.r') # modified to allow adding title

# pairs plot to show correlation
file_name <- paste(out_path, 'last_param_pairs_raw.png', sep="")
png(file_name, height=9, width=9, units='in', res=300)
temp <- pairs(fit, pars=stan_pars_raw, main=setname, 
              xlim=c(0,1), ylim=c(0,1), xbreaks=c(0,0.5,1), ybreaks=c(0,0.5,1))
dev.off()

# pairs plot to show correlation
file_name <- paste(out_path, 'last_param_pairs.png', sep="")
png(file_name, height=9, width=9, units='in', res=300)
temp <- pairs(fit, pars=stan_pars, main=setname)
dev.off()

# pairs plot to show correlation
some_vars <- c('medb0', 'meda1', 'slowb0', 'slowa1', 
               'medBFImax', 'medrec', 'slowBFImax', 'slowrec' )
file_name <- paste(out_path, 'last_param_pairs_paper.png', sep="")
png(file_name, height=9, width=9, units='in', res=300)
temp <- pairs(fit, pars=some_vars, main=setname, priors=priors)
dev.off()

# diagnostic plots to check convergence 
some_vars <- c('medb0', 'meda1', 'slowb0', 'slowa1', 
               'chem1fast0', 'chem1med0', 'chem1slow0',
               'chem2fast0', 'chem2med0', 'chem2slow0',
               'chem1fast1', 'chem1med1', 'chem1slow1',
               'chem2fast1', 'chem2med1', 'chem2slow1')
temp <- rstan::traceplot(fit, pars=some_vars) + labs(title=setname) # crashes in 32-bit R!
# print(temp)
file_name <- paste(out_path, 'last_param_trace.png', sep="")
save_plot(file_name, temp, base_height=9, base_width=12)

some_vars <- c('medb0raw', 'medd1raw', 'slowb0raw', 'slowd1raw',
               'chem1fastraw0', 'chem1medraw0', 'chem1slowraw0',
               'chem2fastraw0', 'chem2medraw0', 'chem2slowraw0',
               'chem1fastraw1', 'chem1medraw1', 'chem1slowraw1',
               'chem2fastraw1', 'chem2medraw1', 'chem2slowraw1')
temp <- traceplot(fit, pars=some_vars) + labs(title=setname)
# print(temp)
file_name <- paste(out_path, 'last_param_trace_raw.png', sep="")
save_plot(file_name, temp, base_height=9, base_width=12)

some_vars <- c('fastflowpc', 'medflowpc', 'slowflowpc',
               'fastTPloadpc', 'medTPloadpc', 'slowTPloadpc', 
               'fastTNloadpc', 'medTNloadpc', 'slowTNloadpc', 
               'nTP', 'nTN')
temp <- traceplot(fit, pars=some_vars) + labs(title=setname)
# print(temp)
file_name <- paste(out_path, 'last_partitioning_trace.png', sep="")
save_plot(file_name, temp, base_height=9, base_width=12)

some_vars <- c('llchem1', 'llchem2', 'vllchem1', 'vllchem2', 
               'grmse1', 'grmse2', 'vgrmse1', 'vgrmse2',
               'auto1', 'auto2', 'vauto1', 'vauto2')
temp <- traceplot(fit, pars=some_vars) + labs(title=setname)
# print(temp)
file_name <- paste(out_path, 'last_fit_trace.png', sep="")
save_plot(file_name, temp, base_height=9, base_width=12)

# use this if necessary to kill processes that didn't terminate
# kill_all_Rscript_s()

# plots
source('box_plots15.r')
source('trace_plots13.r')
