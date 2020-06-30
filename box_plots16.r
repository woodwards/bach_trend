# attempt to make box plots
# box_plots11.r - add shaded bars to shown medians on some plots
# box_plots12.r - add totals to several boxplot, reduce whitespace

library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggthemes)

# path for output
# out_path <- 'run - eckhardt_priors_narrow/'
print(out_path)

# years of simulation
nyears <- 15
print(paste(nyears, "years"))

# functions
"%notin%" <- function(x,y)!("%in%"(x,y))

# read in data
print("reading data")
source('read_data9.r')
nruns <- nrow(runlist)

# read in model samples
print("organising data")
rows <- 1:nruns
datalist <- list()
for (i in rows) {

  # write header
  data_file_name <- paste(runlist$catchfile[i], '_data.dat', sep='')
  opt_file_name <- runlist$optfile[i]

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
  setseq <- aalloptions$setseq
  
  # read samples
  temp <- read_rds(paste(out_path, setname, '_samples.rds', sep=''))
  
  # files are missing for Ot(b) and Wp(b), you can use this to create them
  if (FALSE) {
  temp[, ] <- NA
  temp$setname <- 'Ot(b)'
  write_rds(temp, paste(out_path, 'Ot(b)', '_samples.rds', sep=''))
  temp$setname <- 'Wp(b)'
  write_rds(temp, paste(out_path, 'Wp(b)', '_samples.rds', sep=''))
  }
  
  temp$setseq <- setseq
  datalist[[i]] <- temp # build a list of dataframes

}
sampledata <- bind_rows(datalist) 
sampledata$setseqf <- factor(sampledata$setseq, levels=c(1:nruns))
sampledata$catchname <- substr(sampledata$setname, 1, 2) 
sampledata$catchnamef <- factor(sampledata$catchname)

# define boxplot stats
# https://github.com/tidyverse/ggplot2/issues/898
custombox1 <- function(y) { # for error bar
  data.frame(ymin=quantile(y,0.025),
             lower=quantile(y,0.25),
             middle=quantile(y,0.5),
             upper=quantile(y,0.75),
             ymax=quantile(y,0.975),
             y=y,
             width=0.4)[1, ]
}
custombox2 <- function(y) { # for box
  data.frame(ymin=quantile(y,0.025),
             lower=quantile(y,0.25),
             middle=quantile(y,0.5),
             upper=quantile(y,0.75),
             ymax=quantile(y,0.975),
             y=y,
             width=0.4)[1, ]
}
custombox3 <- function(y) { # for error bar
  data.frame(ymin=quantile(y,0.025),
             lower=quantile(y,0.25),
             middle=quantile(y,0.5),
             upper=quantile(y,0.75),
             ymax=quantile(y,0.975),
             y=y,
             width=0.25)[1, ]
}
custombox4 <- function(y) { # for box
  data.frame(ymin=quantile(y,0.025),
             lower=quantile(y,0.25),
             middle=quantile(y,0.5),
             upper=quantile(y,0.75),
             ymax=quantile(y,0.975),
             y=y,
             width=0.3)[1, ]
}

# test inputs
varname <- 'medb0'
ylabel <- 'Medium b0\n'
ybreaks <- seq(0,1,0.2)
yref <- 0
ytrans <- "identity"
median_fill <- NA
median_colour <- NA
median_alpha <- 1

# handle number of runs
if (nrow(runlist)==24){
  xlabels <- c('','Wh','',
               '','Wp','',
               '','Po','',
               '','Ot','',
               '','Ta','',
               '','Pu','',
               '','Wt','',
               '','Pi','')
  xint <- c(3.5,6.5,9.5,12.5,15.5,18.5,21.5)
  barwidth <- 2.2
  barfilter <- seq(2,nruns,3)
}else {
  xlabels <- c('Wh',
               'Wp',
               'Po',
               'Ot',
               'Ta',
               'Pu',
               'Wt',
               'Pi')
  xint <- c(1.5,2.5,3.5,4.5,5.5,6.5,7.5)
  barwidth <- 0.7
  barfilter <- 1:nruns
}

# box plot. do it using a function
bplot <- function(sampledata, varname, ylabel, ybreaks, yref=0, ytrans="identity"){
  ggplot(data=sampledata) +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox1, geom='errorbar') +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    # https://stackoverflow.com/questions/17479793/changing-bar-width-when-using-stat-summary-with-ggplot
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox2, geom='boxplot') + # width=? is set in custombox
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    scale_x_discrete(labels=xlabels) 
}

# testing
p1 <- bplot(sampledata, 'medb0', 'Medium b0\n', seq(0,1,0.2), c(0,0.5))
print(p1)

# colours
choose <- c(3,5,7,5,3)
tpcol <- brewer.pal(9,"OrRd")[choose]
tncol <- brewer.pal(9,"YlGn")[choose]
fcol <- brewer.pal(9,"YlOrBr")[choose-1]
mcol <- brewer.pal(9,"PuBu")[choose+1]
scol <- brewer.pal(9,"RdPu")[choose+1]

# box plot with medians (which make it a lot slower)
bplot1 <- function(sampledata, varname, ylabel, ybreaks, yref=0, box_fill="white",
                  ytrans="identity", median_fill=NA, median_colour=NA, median_alpha=1){
  # medians <- sampledata[, c(varname, "setname", "setseq", "setseqf")] %>%
  #   rename(value=varname) %>%
  #   group_by(catchnamef) %>%
  #   mutate(median=median(value, na.rm=TRUE))
  if (is.na(median_fill)){
    sampledata2 <- sampledata %>%
      rename(varname=varname) %>%
      select(varname, catchnamef, setseqf, setseq) %>%
      mutate(median=NA_real_)
  } else {
    # this is very slow!!! because lots of copies of median??
    sampledata2 <- sampledata %>%
      rename(varname=varname) %>%
      select(varname, catchnamef, setseqf, setseq) %>%
      group_by(catchnamef) %>%
      mutate(median=median(varname, na.rm=TRUE)) %>%
      arrange(setseq) 
    sampledata2$median[sampledata2$setseq==lag(sampledata2$setseq)] <- NA_real_ # discard dupes
    sampledata2$median[sampledata2$setseq %notin% barfilter ] <- NA_real_ # discard dupes
  }
  ggplot(data=sampledata2) +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_col(mapping=aes(x=setseq, y=median), position="dodge", size=1, width=barwidth, 
             alpha=median_alpha, fill=median_fill, colour=median_colour) +
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y='varname'),
                 fun.data=custombox1, geom='errorbar') +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y='varname'),
                 fun.data=custombox2, geom='boxplot', fill=box_fill) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    # geom_segment(mapping=aes(x=setseq-0.5, xend=setseq+0.5, y=median, yend=median), linetype=2, size=1, colour=median_colour) +
    scale_x_discrete(labels=xlabels, expand=c(0, 0)) 
}

# testing
p1 <- bplot1(sampledata, 'medb0', 'Medium b0\n', seq(0,1,0.2), c(0,0.5), median_fill=mcol[3], median_colour=mcol[3], median_alpha=0.2)
print(p1)

# pars
print("making pars box")
p1 <- bplot(sampledata, 'medb0', 'Medium b0\n', seq(0,1,0.2), c(0,0.9))
p2 <- bplot(sampledata, 'slowb0', 'Slow b0\n', seq(0,0.005,0.001), c(0,0.01))
p3 <- bplot(sampledata, 'meda1', 'Medium a1\n', seq(0,1,0.2), c(0.1,0.99))
p4 <- bplot(sampledata, 'slowa1', 'Slow a1\n', seq(0.995,1,0.001), c(0.99,0.9999))
plotbox <- plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'pars', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=8)

# Eckhardt pars
print("making eckpars box")
p1 <- bplot(sampledata, 'medBFImax', 'Medium BFImax\n', seq(0,1,0.2), c(0,1))
p2 <- bplot(sampledata, 'slowBFImax', 'Slow BFImax\n', seq(0,1,0.2), c(0,1))
p3 <- bplot(sampledata, 'medrec', 'Medium k\n', seq(0.9,1,0.02), c(0.5,1))
p4 <- bplot(sampledata, 'slowrec', 'Slow k\n', seq(0.999,1,0.0002), c(0.99,1))
plotbox <- plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'eckpars', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=8)

# # synthetic pars (trying to reduce unc)
# print("making synth box")
# sampledata$synth_m1 <- log(1-sampledata$meda1)
# sampledata$synth_s1 <- log(1-sampledata$slowa1)
# sampledata$synth_m2 <- sampledata$medb0 * sampledata$meda1
# sampledata$synth_s2 <- sampledata$slowb0 * sampledata$slowa1
# p1 <- bplot(sampledata, 'synth_m1', 'Medium Synth 1\n', seq(-6,0,1), c(0,1))
# p2 <- bplot(sampledata, 'synth_s1', 'Slow Synth 1\n', seq(-10,-5,1), c(0,1))
# p3 <- bplot(sampledata, 'synth_m2', 'Medium Synth 2\n', seq(0,1,0.1), c(0,1))
# p4 <- bplot(sampledata, 'synth_s2', 'Slow Synth 2\n', seq(0,0.005,0.001), c(0,1))
# plotbox <- plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
# # print(plotbox)
# file_name <- paste(out_path, 'box_', 'synth', '.png', sep="")
# save_plot(file_name, plotbox, base_height=6, base_width=8)

# concs
print("making conc0 box")
p1 <- bplot1(sampledata, 'chem1fast0', expression('Fast TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'chem1med0', expression('Medium TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'chem1slow0', expression('Slow TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'chem2fast0', expression('Fast TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'chem2med0', expression('Medium TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'chem2slow0', expression('Slow TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'conc0', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

print("making conc1 box")
p1 <- bplot1(sampledata, 'chem1fast1', expression('Fast TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'chem1med1', expression('Medium TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'chem1slow1', expression('Slow TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'chem2fast1', expression('Fast TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'chem2med1', expression('Medium TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'chem2slow1', expression('Slow TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'conc1', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

print("making concTP box")
p1 <- bplot1(sampledata, 'chem1fast0', expression('Fast TP Conc0'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'chem1med0', expression('Medium TP Conc0'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'chem1slow0', expression('Slow TP Conc0'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'chem1fast1', expression('Fast TP Conc1'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'chem1med1', expression('Medium TP Conc1'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'chem1slow1', expression('Slow TP Conc1'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'concTP', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

print("making concTN box")
p1 <- bplot1(sampledata, 'chem2fast0', expression('Fast TN Conc0'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'chem2med0', expression('Medium TN Conc0'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'chem2slow0', expression('Slow TN Conc0'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'chem2fast1', expression('Fast TN Conc1'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'chem2med1', expression('Medium TN Conc1'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'chem2slow1', expression('Slow TN Conc1'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'concTN', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

# side by side two variables
bplot2 <- function(sampledata, varname0, varname1, ylabel, ybreaks, yref=0, box_fill="white",
                   ytrans="identity", median_fill=NA, median_colour=NA, median_alpha=1){
  if (is.na(median_fill)){
    sampledata2 <- sampledata %>%
      rename(varname0=varname0, varname1=varname1) %>%
      select(varname0, varname1, catchnamef, setseqf, setseq) %>%
      mutate(median0=NA_real_, median1=NA_real_)
  } else {
    # this is very slow!!! because lots of copies of median??
    sampledata2 <- sampledata %>%
      rename(varname0=varname0, varname1=varname1) %>%
      select(varname0, varname1, catchnamef, setseqf, setseq) %>%
      group_by(catchnamef) %>%
      mutate(median0=median(varname0, na.rm=TRUE), median1=median(varname1, na.rm=TRUE)) %>%
      arrange(setseq) 
    sampledata2$median0[sampledata2$setseq==lag(sampledata2$setseq)] <- NA_real_ # discard dupes
    sampledata2$median0[sampledata2$setseq %notin% barfilter ] <- NA_real_ # discard dupes
    sampledata2$median1[sampledata2$setseq==lag(sampledata2$setseq)] <- NA_real_ # discard dupes
    sampledata2$median1[sampledata2$setseq %notin% barfilter ] <- NA_real_ # discard dupes
  }
  ggplot(data=sampledata2) +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_col(mapping=aes(x=setseq-0.25, y=median0), position="dodge", size=1, width=barwidth/2, 
             alpha=median_alpha, fill=median_fill, colour=median_colour) +
    geom_col(mapping=aes(x=setseq+0.25, y=median1), position="dodge", size=1, width=barwidth/2, 
             alpha=median_alpha, fill=median_fill, colour=median_colour) +
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(mapping=aes(group=setseq, x=setseq-0.25, y=varname0), fun.data=custombox3, geom='errorbar') +
    stat_summary(mapping=aes(group=setseq, x=setseq-0.25, y=varname0), fun.data=custombox4, geom='boxplot', fill=box_fill) +
    stat_summary(mapping=aes(group=setseq, x=setseq+0.25, y=varname1), fun.data=custombox3, geom='errorbar') +
    stat_summary(mapping=aes(group=setseq, x=setseq+0.25, y=varname1), fun.data=custombox4, geom='boxplot', fill=box_fill) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    # geom_segment(mapping=aes(x=setseq-0.5, xend=setseq+0.5, y=median, yend=median), linetype=2, size=1, colour=median_colour) +
    scale_x_continuous(breaks=seq_along(xlabels), labels=xlabels, expand=c(0, 0)) 
}

print("making dconc box")
p1 <- bplot2(sampledata, 'chem1fast0', 'chem1fast1', expression('Fast TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot2(sampledata, 'chem1med0', 'chem1med1', expression('Medium TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot2(sampledata, 'chem1slow0', 'chem1slow1', expression('Slow TP Conc'~(mg~L^{-1})*'\n'), seq(0,0.8,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot2(sampledata, 'chem2fast0', 'chem2fast1', expression('Fast TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot2(sampledata, 'chem2med0', 'chem2med1', expression('Medium TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot2(sampledata, 'chem2slow0', 'chem2slow1', expression('Slow TN Conc'~(mg~L^{-1})*'\n'), seq(0,6,1), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'dconc', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

# loads
sampledata <- sampledata %>%
  mutate(
    fastflowmm = fastflow/area*86.4*365.25,
    medflowmm = medflow/area*86.4*365.25,
    slowflowmm = slowflow/area*86.4*365.25,
    totalflowmm = fastflowmm+medflowmm+slowflowmm,
    fastTPloadkg = fastTPload/area*10/nyears, # convert from t to kg/ha/y
    medTPloadkg = medTPload/area*10/nyears,
    slowTPloadkg = slowTPload/area*10/nyears,
    totalTPloadkg = fastTPloadkg+medTPloadkg+slowTPloadkg,
    fastTNloadkg = fastTNload/area*10/nyears,
    medTNloadkg = medTNload/area*10/nyears,
    slowTNloadkg = slowTNload/area*10/nyears,
    totalTNloadkg = fastTNloadkg+medTNloadkg+slowTNloadkg,
    fastTPloadpc = fastTPloadkg/(fastTPloadkg+medTPloadkg+slowTPloadkg),
    medTPloadpc = medTPloadkg/(fastTPloadkg+medTPloadkg+slowTPloadkg),
    slowTPloadpc = slowTPloadkg/(fastTPloadkg+medTPloadkg+slowTPloadkg),
    fastTNloadpc = fastTNloadkg/(fastTNloadkg+medTNloadkg+slowTNloadkg),
    medTNloadpc = medTNloadkg/(fastTNloadkg+medTNloadkg+slowTNloadkg),
    slowTNloadpc = slowTNloadkg/(fastTNloadkg+medTNloadkg+slowTNloadkg)
  )
#
print("making totals boxes")
tcol <- tncol[3]
varname <- 'totalTPloadkg'
pt1 <- bplot1(sampledata, 'totalTPloadkg', expression('Total TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1.6,0.2), median_colour=NA, median_fill=tcol, median_alpha=0.2)
varname <- 'totalTNloadkg'
pt2 <- bplot1(sampledata, 'totalTNloadkg', expression('Total TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=tcol, median_alpha=0.2)
varname <- 'totalflowmm'
pt3 <- bplot1(sampledata, 'totalflowmm', expression('Total Flow'~(mm~y^{-1})~'\n'), seq(0,1200,200), median_colour=NA, median_fill=tcol, median_alpha=0.2)
plotbox <- pt3
file_name <- paste(out_path, 'box_', 'totalflow', '.png', sep="")
save_plot(file_name, plotbox, base_height=3, base_width=4)
plotbox <- plot_grid(pt1, pt2, nrow=2, align="hv")
file_name <- paste(out_path, 'box_', 'totalload', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=4)
plotbox <- plot_grid(pt3, pt1, pt2, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'totalall', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=4)
# plotbox <- plot_grid(p1, p2, nrow=2, align="hv")
# file_name <- paste(out_path, 'box_', 'totalload', '.png', sep="")
# save_plot(file_name, plotbox, base_height=6, base_width=4)
# file_name <- paste(out_path, 'box_', 'totalflow', '.png', sep="")
# save_plot(file_name, p3, base_height=3, base_width=4)

print("making TPload box")
p1 <- bplot1(sampledata, 'fastTPloadkg', expression('Fast TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'medTPloadkg', expression('Medium TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'slowTPloadkg', expression('Slow TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'fastTPloadpc', 'Fast TP Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'medTPloadpc', 'Medium TP Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
varname <- 'slowTPloadpc'
p6 <- bplot1(sampledata, 'slowTPloadpc', 'Slow TP Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'TPload', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)
# plotbox <- plot_grid(pt1, NULL, p1, p4, p2, p5, p3, p6, nrow=4, align="hv")
# file_name <- paste(out_path, 'box_', 'TPload', '.png', sep="")
# save_plot(file_name, plotbox, base_height=12, base_width=8)

print("making TNload box")
p1 <- bplot1(sampledata, 'fastTNloadkg', expression('Fast TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'medTNloadkg', expression('Medium TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'slowTNloadkg', expression('Slow TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'fastTNloadpc', 'Fast TN Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'medTNloadpc', 'Medium TN Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'slowTNloadpc', 'Slow TN Yield Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'TNload', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)
# plotbox <- plot_grid(pt2, NULL, p1, p4, p2, p5, p3, p6, nrow=4, align="hv")
# file_name <- paste(out_path, 'box_', 'TNload', '.png', sep="")
# save_plot(file_name, plotbox, base_height=12, base_width=8)

print("making BOTHload box")
p1 <- bplot1(sampledata, 'fastTPloadkg', expression('Fast TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'medTPloadkg', expression('Medium TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'slowTPloadkg', expression('Slow TP Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
p4 <- bplot1(sampledata, 'fastTNloadkg', expression('Fast TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'medTNloadkg', expression('Medium TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'slowTNloadkg', expression('Slow TN Yield '~(kg~ha^{-1}~y^{-1})*'\n'), seq(0,20,5), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'BOTHload', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)
# plotbox <- plot_grid(pt2, NULL, p1, p4, p2, p5, p3, p6, nrow=4, align="hv")
# file_name <- paste(out_path, 'box_', 'TNload', '.png', sep="")
# save_plot(file_name, plotbox, base_height=12, base_width=8)

# flowpc
print("making flow box")
p4 <- bplot1(sampledata, 'fastflowpc', 'Fast Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'medflowpc', 'Medium Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'slowflowpc', 'Slow Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
# plotbox <- plot_grid(p1, p2, p3, nrow=3, align="hv")
# file_name <- paste(out_path, 'box_', 'flowpc', '.png', sep="")
# save_plot(file_name, plotbox, base_height=9, base_width=4)
# flow
p1 <- bplot1(sampledata, 'fastflowmm', expression('Fast Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'medflowmm', expression('Medium Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'slowflowmm', expression('Slow Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
# plotbox <- plot_grid(p1, p2, p3, nrow=3, align="hv")
# file_name <- paste(out_path, 'box_', 'flowmm', '.png', sep="")
# save_plot(file_name, plotbox, base_height=9, base_width=4)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'flow', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)
# plotbox <- plot_grid(pt3, NULL, p1, p4, p2, p5, p3, p6, nrow=4, align="hv")
# file_name <- paste(out_path, 'box_', 'flow', '.png', sep="")
# save_plot(file_name, plotbox, base_height=12, base_width=8)

# this version allows different ref lines for auto
bplot2 <- function(sampledata, varname, ylabel, ybreaks, yref1, yref2){
  # medfn <- function(x){
  #   med <- median(x)
  #   if (med==0) med=NA
  #   return(med)
  # }
  ggplot() +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_hline(yintercept=0, linetype=1, colour='black') +
    geom_path(aes(x=1:24, y=yref1), linetype=2, colour='black') +
    geom_path(aes(x=1:24, y=yref2), linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(data=sampledata, mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox1, geom='errorbar') +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    stat_summary(data=sampledata, mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox2, geom='boxplot') +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0)) +
    scale_x_discrete(labels=xlabels) 
}
# auto
print("making auto box FIXME bug!")
yscale <- c(-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8)
refline1 <- rep(-0.25303,24)
refline2 <- rep(0.25303,24)
refline1[c(4,10)] <- -0.331294
refline2[c(4,10)] <- 0.331294
p1 <- bplot2(sampledata, 'auto1', 'R(1) TP Calib\n', yscale, refline1, refline2)
# print(p1)
refline1 <- rep(-0.25303,24)
refline2 <- rep(0.25303,24)
p2 <- bplot2(sampledata, 'vauto1', 'R(1) TP Valid\n', yscale, refline1, refline2)
p3 <- bplot2(sampledata, 'auto2', 'R(1) TN Calib\n', yscale, refline1, refline2)
p4 <- bplot2(sampledata, 'vauto2', 'R(1) TN Valid\n', yscale, refline1, refline2)
plotbox <- plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'auto', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=8)

# gmrse 
print("making grmse box")
sampledata$vgrmse1[sampledata$vgrmse1==0] <- NA
sampledata$vgrmse2[sampledata$vgrmse2==0] <- NA
p1 <- bplot(sampledata, 'grmse1', expression('GRMSE TP Calib'~(mg~L^{-1})*'\n'), seq(0,0.07,0.01), 0.02)
p2 <- bplot(sampledata, 'vgrmse1', expression('GRMSE TP Valid'~(mg~L^{-1})*'\n'), seq(0,0.07,0.01), 0.02)
p3 <- bplot(sampledata, 'grmse2', expression('GRMSE TN Calib'~(mg~L^{-1})*'\n'), seq(0,0.7,0.1), 0.2)
p4 <- bplot(sampledata, 'vgrmse2', expression('GRMSE TN Valid'~(mg~L^{-1})*'\n'), seq(0,0.7,0.1), 0.2)
# plotbox <- plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
plotbox <- plot_grid(p1, p3, nrow=2, align="hv")
# print(plotbox)
file_name <- paste(out_path, 'box_', 'grmse', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=4)

# bar plot function
barplot <- function(sampledata, varname, ylabel, ybreaks, yref1=NA, yref2=NA, ytrans="identity"){
  ggplot() +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_bar(data=sampledata, mapping=aes_string(x='setseqf', y=varname), stat="identity", fill="cornflowerblue", width=0.7) +
    geom_hline(yintercept=yref1, linetype=2, colour='black') +
    geom_hline(yintercept=yref2, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    #	  stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, expand=c(0, 0), trans=ytrans) +
    coord_cartesian(ylim=c(min(ybreaks), max(ybreaks))) + # allows bars to go off the page
    scale_x_discrete(labels=xlabels) 
}

# gelman r 
print("making gelman box")
runfile <- paste(out_path,'runrecord.tsv', sep='')
runrecord <- read_tsv(runfile, col_types=cols())
runrecord$setseqf <- factor(runrecord$setseq, levels=c(1:nruns))
p1 <- barplot(runrecord, 'gelmanr', 'Gelman R\n', seq(0,max(runrecord$gelmanr)+0.5,0.5), 1.1, 1.0)
temp <- ceiling(max(runrecord$elapsed, na.rm=TRUE))
p2 <- barplot(runrecord, 'elapsed', 'Elapsed (h)\n', seq(0,temp,0.5), 0, 0)
plotbox <- plot_grid(p1, p2, nrow=2, align="v")
print(plotbox)
file_name <- paste(out_path, 'box_', 'gelmanr', '.png', sep="")
save_plot(file_name, plotbox, base_height=6, base_width=4)

# redefine boxplot stats
custombox <- function(y) {
  data.frame(ymin=quantile(y,0.0),
             lower=quantile(y,0.25),
             middle=quantile(y,0.5),
             upper=quantile(y,0.75),
             ymax=quantile(y,1.0))
}

# annual flow
bplot3 <- function(sampledata, varname, ylabel, ybreaks, yref, ytrans="identity"){
  ggplot() +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(data=sampledata, mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox1, geom='errorbar') +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    stat_summary(data=sampledata, mapping=aes_string(group='setseqf', x='setseqf', y=varname),
                 fun.data=custombox2, geom='boxplot') +
    geom_point(data=sampledata, mapping=aes_string(x='setseqf', y=varname), size=1) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname),
    #              fun.y=mean, geom='point', pch=21, size=1, fill="white", colour="black") +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    scale_x_discrete(labels=xlabels) 
}

# annual flow alternative version
bplot4 <- function(sampledata, varname, ylabel, ybreaks, yref, ytrans="identity"){
  ggplot() +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname),
    #              fun.data=custombox, geom='errorbar', width=0.35) +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname),
    #              fun.data=custombox, geom='boxplot') +
    geom_point(data=sampledata, mapping=aes_string(x='setseqf', y=varname), size=1) +
    geom_point(data=sampledata, mapping=aes_string(x='setseqf', y=varname), size=2, pch=1, stat="summary", fun.y="mean") +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname),
    #              fun.y=mean, geom='point', pch=21, size=1, fill="white", colour="black") +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    scale_x_discrete(labels=xlabels) 
}

# annual flow
print("making annual_flow box")
data2 <- data %>%
  filter(period!="None") %>%
  group_by(setname, year, setseq) %>%
  summarize(annualmm=sum(flow2),
            meancumec=mean(flow),
            meanTP=mean(TP, na.rm=TRUE),
            meanTN=mean(TN, na.rm=TRUE))
data2$setseqf <- factor(data2$setseq, levels=c(1:nruns))
p1 <- bplot4(data2, 'annualmm', 'Annual Flow', seq(0,1400,200))
print(p1)
p2 <- bplot4(data2, 'meanTP', 'Annual Mean TP', seq(0,0.3,0.1))
print(p2)
p3 <- bplot4(data2, 'meanTN', 'Annual Mean TN', seq(0,3,1))
print(p3)
plotbox <- plot_grid(p2, p3, p1, nrow=3, align="v")
print(plotbox)
file_name <- paste(out_path, 'box_', 'annual_flow', '.png', sep="")
save_plot(file_name, p1, base_height=3, base_width=4)
file_name <- paste(out_path, 'box_', 'annual_summary', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=4)

data3 <- data %>%
  filter(period!="None") %>%
  group_by(shortname) %>%
  summarize(avmm=sum(flow2)/length(unique(year)),
            avcumec=mean(flow))

# histograms
# greens 2 # https://www.color-hex.com/color-palette/5016
xpale <- "#F0F7DA"
xlight <- "#C9DF8A"
xmid <- "#77AB59"
xdark <- "#36802D"
xdata <- "#043927" # https://graf1x.com/shades-of-green-color-palette-html-hex-rgb-code/
xaxis <- "#999999"
xgrey <- "grey"
prior_df <- vector("list", length(priortab$parname))
i <- 1
for (i in seq_along(priortab$parname)) {
  key <- priortab$parname[i]
  x <- seq(priortab$parmin[i], priortab$parmax[i], length.out = 101) * priortab$parscale[i]
  y <- dnorm(x, priortab$parmean[i] * priortab$parscale[i], priortab$parsd[i] * priortab$parscale[i])
  x <- case_when(
    key == "medd1raw" ~ x + (-0.1*log(1-0.1)) * 10,
    key == "slowd1raw" ~ x + (-0.1*log(1-0.99)) * 10,
    TRUE ~ x)
  prior_df[[i]] <- tibble(key = factor(key, levels = priortab$parname), x = x, y = y)
}
prior_df <- bind_rows(prior_df) 
my_pretty_breaks <- function(n = 5, ...) {
  n_default <- n
  function(x, n = n_default) {
    minx <- min(x)
    maxx <- max(x)
    midx <- (minx + maxx) / 2
    x2 <- midx + (x - midx) * 0.9
    breaks <- pretty(x2, n, ...)
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
}
rows <- 1:nruns
for (i in rows) {
  runlisti <- runlist %>% 
    filter(setseq == i)
  print(paste("making par histogram", runlisti$setname))
  sampledatai <- sampledata %>% 
    filter(setname == runlisti$setname) %>% 
    select_at(priortab$parname)
  scalelist <- setNames(priortab$parscale, priortab$parname)
  post_df <- sampledatai %>% 
    pivot_longer(everything(), names_to = "key", values_to = "value") %>% 
    mutate(
      key = factor(key, levels = priortab$parname),
      value = value * scalelist[key],
      value = case_when(
        key == "medd1raw" ~ value + (-0.1*log(1-0.5)) * 10,
        key == "slowd1raw" ~ value + (-0.1*log(1-0.99)) * 10,
        TRUE ~ value)
      ) 
  labellist <- setNames(priortab$parlabel, priortab$parname)
  getlabels <- function(x){
    unname(labellist[x])
  }
  plot1 <- ggplot(data = prior_df) +
    labs(title = runlisti$catchname, x = "", y = "") +
    geom_line(mapping = aes(x = x, y = y), colour = xmid, size = 1) +
    geom_histogram(data = post_df, mapping = aes(x = value, y = ..density..), fill = xlight, bins = 30) +
    geom_line(mapping = aes(x = x, y = y), colour = xmid, size = 1) +
    theme_few() +
    theme(panel.spacing.x = unit(9, "mm"), plot.margin = unit(c(0, 5, 0, 0), "mm")) +
    # theme(axis.text=element_text(size=9)) +
    scale_x_continuous(expand = expand_scale(0, 0), breaks = my_pretty_breaks(n = 2, min.n = 2)) +
    scale_y_continuous(expand = expand_scale(0, 0), breaks = NULL) +
    facet_wrap(vars(key), nrow = 4, scales = "free", labeller = labeller(key = getlabels))
  # print(plot1)
  png(paste0(out_path, "/", runlisti$setname, "_histogram.png"),
      width = 210 * 1.5, height = 210 * 1, units = "mm",
      type = "windows", res = 300 # FIXME res = 600
  )
  print(plot1)
  dev.off()
}
  
# get old sampledata for comparison
print("reading old sampledata")
oldsampledata <- readRDS(paste0(out_path,"oldsampledata.rds"))

custombox1top <- function(y) { # for error bar
  data.frame(ymin=quantile(y,0.025),
             lower=quantile(y,0.25),
             middle=quantile(y,0.25),
             upper=quantile(y,0.25),
             ymax=quantile(y,0.25),
             y=y,
             width=0.4)[1, ]
}
custombox1bot <- function(y) { # for error bar
  data.frame(ymin=quantile(y,0.75),
             lower=quantile(y,0.75),
             middle=quantile(y,0.75),
             upper=quantile(y,0.75),
             ymax=quantile(y,0.975),
             y=y,
             width=0.4)[1, ]
}

# box plot with medians (which make it a lot slower)
bplot1 <- function(sampledata, varname, ylabel, ybreaks, yref=0, box_fill=NA,
                   ytrans="identity", median_fill=NA, median_colour=NA, median_alpha=1){
  # medians <- sampledata[, c(varname, "setname", "setseq", "setseqf")] %>%
  #   rename(value=varname) %>%
  #   group_by(catchnamef) %>%
  #   mutate(median=median(value, na.rm=TRUE))
  if (is.na(median_fill)){
    sampledata2 <- sampledata %>%
      rename(varname=varname) %>%
      select(varname, catchnamef, setseqf, setseq) %>%
      mutate(median=NA_real_)
  } else {
    # this is very slow!!! because lots of copies of median??
    sampledata2 <- sampledata %>%
      rename(varname=varname) %>%
      select(varname, catchnamef, setseqf, setseq) %>%
      group_by(catchnamef) %>%
      mutate(median=median(varname, na.rm=TRUE)) %>%
      arrange(setseq) 
    sampledata2$median[sampledata2$setseq==lag(sampledata2$setseq)] <- NA_real_ # discard dupes
    sampledata2$median[sampledata2$setseq %notin% barfilter ] <- NA_real_ # discard dupes
    oldsampledata2 <- oldsampledata %>%
      rename(varname=varname) %>%
      select(varname, catchnamef, setseqf, setseq) %>%
      group_by(catchnamef) %>%
      mutate(median=median(varname, na.rm=TRUE)) %>%
      arrange(setseq) %>% 
      mutate(
        # setseq = (setseq - 1) %/% 3 + 1,
        setseqf = factor(setseq)
        )
    oldsampledata2$median[oldsampledata2$setseq==lag(oldsampledata2$setseq)] <- NA_real_ # discard dupes
    oldsampledata2$median[oldsampledata2$setseq %notin% barfilter ] <- NA_real_ # discard dupes
  }
  ggplot(data=sampledata2) +
    labs(title='', y=ylabel, x='') +
    theme_cowplot(font_size=10) +
    theme(axis.ticks.x=element_blank(), 
          # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.title=element_blank()) +
    panel_border(colour='black') +  
    # geom_col(data = oldsampledata2,
    #          mapping=aes(x=setseq, y=median), position="dodge", size=1, width=barwidth,
    #          alpha=median_alpha, fill=median_fill, colour=median_colour) +
    geom_violin(data = oldsampledata2,
             mapping=aes(x=(setseq+1)/3, y=varname, group=setseqf), width=0.8,
             alpha=0.4, fill=median_fill, colour=median_colour) +
    geom_hline(yintercept=yref, linetype=2, colour='black') +
    geom_vline(xintercept=xint, colour='grey') +
    # stat_boxplot(data=sampledata, mapping=aes_string(x='setseqf', y=varname), 
    #             geom="errorbar", width=0.35) +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y='varname'),
                 fun.data=custombox1top, geom='errorbar') +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y='varname'),
                 fun.data=custombox1bot, geom='errorbar') +
    # geom_boxplot(outlier.size=0.5, notch=FALSE, outlier.shape=NA) +
    stat_summary(mapping=aes_string(group='setseqf', x='setseqf', y='varname'),
                 fun.data=custombox2, geom='boxplot', fill=box_fill) +
    # stat_summary(data=sampledata, mapping=aes_string(x='setseqf', y=varname), fun.y=median, geom='point', pch=3, size=1) +
    scale_y_continuous(breaks=ybreaks, limits=c(min(ybreaks), max(ybreaks)), expand=c(0, 0), trans=ytrans) +
    # geom_segment(mapping=aes(x=setseq-0.5, xend=setseq+0.5, y=median, yend=median), linetype=2, size=1, colour=median_colour) +
    scale_x_discrete(labels=xlabels, expand=c(0, 0)) 
}

# flowpc
print("making flow box vs old data")
p4 <- bplot1(sampledata, 'fastflowpc', 'Fast Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p5 <- bplot1(sampledata, 'medflowpc', 'Medium Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p6 <- bplot1(sampledata, 'slowflowpc', 'Slow Flow Fraction\n', seq(0,1,0.2), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
# plotbox <- plot_grid(p1, p2, p3, nrow=3, align="hv")
# file_name <- paste(out_path, 'box_', 'flowpc', '.png', sep="")
# save_plot(file_name, plotbox, base_height=9, base_width=4)
# flow
p1 <- bplot1(sampledata, 'fastflowmm', expression('Fast Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=fcol[3], median_alpha=0.2)
p2 <- bplot1(sampledata, 'medflowmm', expression('Medium Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=mcol[3], median_alpha=0.2)
p3 <- bplot1(sampledata, 'slowflowmm', expression('Slow Flow'~(mm~y^{-1})~'\n'), seq(0,1000,200), median_colour=NA, median_fill=scol[3], median_alpha=0.2)
# plotbox <- plot_grid(p1, p2, p3, nrow=3, align="hv")
# file_name <- paste(out_path, 'box_', 'flowmm', '.png', sep="")
# save_plot(file_name, plotbox, base_height=9, base_width=4)
plotbox <- plot_grid(p1, p4, p2, p5, p3, p6, nrow=3, align="hv")
file_name <- paste(out_path, 'box_', 'flow_comparison', '.png', sep="")
save_plot(file_name, plotbox, base_height=9, base_width=8)

