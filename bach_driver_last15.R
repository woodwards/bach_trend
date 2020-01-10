#
print("Generating diagnostics and charts for last run.")

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
  facet_wrap("key")
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
  facet_wrap("key", scales="free")
print(p2)
dev.off()

# use modified code
source('pairs.stanfit.r') # modified to allow adding title

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


