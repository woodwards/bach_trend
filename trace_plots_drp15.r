# trace plots for bach
# add TF plot
# trace_plots11.r - modify colours to match boxplots

# library(MASS)
library(robustbase)
library(tidyverse)
library(cowplot)
library(scales) # for dates
library(rkt) # Mann-Kendall test
library(lubridate)
library(lfstat)
library(FlowScreen)
library(RColorBrewer)
library(naniar)
library(EGRET) # UGSG trend software

# path for output
# out_path <- 'run - eckhardt_priors_narrow/'
print(out_path)
TPname <- "DRP"
TNname <- "NO3-N"

# choose rows to run
rows <- 1:nruns # to run all sites
ah <- vector("list", nruns) # all hydrographs
# rows <- c(20,21) # to run a subset of sites, useful for testing
# i <- rows[[10]]

#### loop through data sets ####
i <- 2
for (i in rows) {

  # set axis options for each set
  TPbreaks <- seq(0, 0.4, 0.1)
  TPlimits <- range(TPbreaks)
  TNbreaks <- seq(0, 3, 1.0)
  TNlimits <- range(TNbreaks)
  flowbreaks <- seq(0, 3, 1)
  flowlimits <- range(flowbreaks)
  
  dTPlimits <- c(-0.1, 0.1)
  dTPbreaks <- seq(-0.1, 0.1, 0.05)
  dTNlimits <- c(-1, 1)
  dTNbreaks <- seq(-1, 1, 0.5)
  
  # rough estimate for total flow (TF)
  TFlimits <- flowlimits * 1
  TFbreaks <- flowbreaks * 1

  data_file_name <- runlist$catchfile[i]
  opt_file_name <- runlist$datafile[i]
  print(paste(data_file_name, "+", opt_file_name))

  # data
  arun <- runlist[i, ]
  temp1 <- flowdata %>% 
    filter(catch == data_file_name) %>% 
    select(date, flow) 
  temp2 <- concdata %>% 
    filter(catch == data_file_name) %>% 
    mutate(date = as.Date(date)) %>% 
    select(date, all_of(c(arun$chem1i, arun$chem2i)))
  adata <- left_join(temp1, temp2, by = "date")
  
  # assemble run options/data into vectors (?)
  nyears <- arun$nyears
  catchname <- arun$catchname
  setname <- arun$setname
  ii <- arun$startrun:arun$endcalib
  jj <- arun$startcalib:arun$endcalib - arun$startrun + 1
  calibxrange <- jj
  refdate <- arun$date1
  
  date <- adata$date[ii]
  ndate <- length(date)
  flow <- as.vector(adata$flow[ii] / 1000)
  TP <- as.vector(adata[[arun$chem1i]][ii] / 1000)
  TN <- as.vector(adata[[arun$chem2i]][ii] / 1000)
  meanTPcalib <- mean(TP[jj], na.rm=TRUE)
  meanTNcalib <- mean(TN[jj], na.rm=TRUE)
  meanTFcalib <- mean(flow[jj], na.rm=TRUE)
  TP[is.na(TP)] <- -1 # stan doesn't understand NA
  TN[is.na(TN)] <- -1 # stan doesn't understand NA
  
  # add flow info
  adata <- adata %>%
    mutate(
      TF = flow/1000,
      shortname = arun$shortname
    )

  # define date range
  fromdate <- adata$date[arun$startcalib]
  todate <- adata$date[arun$endcalib]
  daterange <- c(fromdate, todate)

  # filter on date
  adata <- adata %>%
    slice(ii) %>% 
    mutate(
      x = 1:n(),
      xdate = date
    )
  
  #### do hydrograph ####

  # get hydrograph
  y0 <- adata %>%
    mutate(pc = "data") %>%
    select(pc, TF, x, xdate) %>%
    filter(x %in% calibxrange) %>%
    mutate(
      flow = TF,
      day = day(xdate),
      month = month(xdate),
      year = year(xdate)
    )
  ldata <- createlfobj(y0, hyearstart = 1) # does baseflow separation
  y0$baseflow <- ldata$baseflow
  bfi <- BFI(ldata)
  print(paste("BFI Fixed    =", bfi))
  y0$baseflow2 <- bf_eckhardt(y0$flow, 0.98, 0.8)
  bfi2 <- sum(y0$baseflow2) / sum(y0$flow)
  print(paste("BFI Eckhardt =", bfi2))

  p0 <- ggplot() +
    labs(title = "", y = "Flow " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
    theme_cowplot() +
    theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
    panel_border(colour = "black") +
    # scale_colour_manual(values=c('orchid','darkorchid','darkmagenta','darkorchid','orchid')) +
    scale_colour_manual(values = c("black", "black", "black", "black", "black")) +
    scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
    scale_y_continuous(expand = c(0, 0), breaks = TFbreaks) +
    coord_cartesian(ylim = TFlimits) +
    # geom_ribbon(data=y0, mapping=aes(x=xdate, ymin=TF-0.2*tci*TF, ymax=TF+0.2*tci*TF), colour='lightgrey', fill='lightgrey') +
    geom_line(data = y0, mapping = aes(x = xdate, y = TF, colour = pc))
  aah <- p0 +
    annotate(geom = "label", x = min(y0$xdate) + 120, y = max(TFlimits) * 0.8, label = setname, size = 4) +
    # geom_line(data=y0, mapping=aes(x=xdate, y=baseflow2), colour="red") +
    labs(y = "")
  # print(aah)

  # now do stuff that only applies to data sets that have chem data
  file_name <- paste(out_path, setname, "_traces.rds", sep = "")
  if (file.exists(file_name)) {

    # read traces
    traces <- read_rds(file_name)
    traces[4:17,] <- NA # conc model not valid here
    traces[4:17+17,] <- NA # conc model not valid here
    
    # write to file
    varlist <- c("flowtotal", "flownonfast", "flowslow",
                 "drpmodel", "no3nmodel", 
                 "drpfast", "drpmed", "drpslow", "no3nfast", "no3nmed", "no3nslow",
                 "ddrpfast", "ddrpmed", "ddrpslow", "dno3nfast", "dno3nmed", "dno3nslow")
    temp <- traces %>% 
      rownames_to_column(var = "rowname") %>% 
      mutate(
        timestep = str_extract(rowname, "(?<=\\[)[0-9]+(?=\\,)"),
        vari = str_extract(rowname, "(?<=\\,)[0-9]+(?=\\])"),
        varname = varlist[as.integer(vari)]
      )
    file_name <- paste(out_path, setname, "_traces.csv", sep = "")
    write_csv(temp, file_name)
    
    # this is probably a better approach, use a single data frame for all variables
    # temp <- traces %>%
    #   select('2.5%', '25%', '50%', '75%', '97.5%') %>%
    #   mutate(name=rownames(traces)) %>%
    #   gather('2.5%', '25%', '50%', '75%', '97.5%', key='pc', value='val') %>%
    #   mutate(var=rep(c('Fast', 'Medium', 'Slow', 'TP', 'TN'), times=ndate*5),
    #          date=rep(calibxrange, each=5, times=5))

    sort <- c(1, 5, 2, 4, 3)
    # tpcol <- c('chartreuse','green','seagreen','green','chartreuse')
    # tncol <- c('pink','red','darkred','red','pink')
    # fcol <- c('orange','firebrick','darkorange','firebrick','orange')
    # mcol <- c('deepskyblue','blue','darkblue','blue','deepskyblue')
    # scol <- c('orchid','darkorchid','darkmagenta','darkorchid','orchid')
    choose <- c(3, 5, 7, 5, 3)
    tpcol <- brewer.pal(9, "OrRd")[choose]
    tncol <- brewer.pal(9, "YlGn")[choose]
    fcol <- brewer.pal(9, "YlOrBr")[choose - 1]
    mcol <- brewer.pal(9, "PuBu")[choose + 1]
    scol <- brewer.pal(9, "RdPu")[choose + 1]

    keeps <- grep(",3]", rownames(traces), invert = FALSE)
    y3 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "Slow") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)

    keeps <- grep(",2]", rownames(traces), invert = FALSE)
    y2 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "Medium") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange) %>%
      mutate(Medium = Medium - y3$Slow)

    keeps <- grep(",1]", rownames(traces), invert = FALSE)
    y1 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "Fast") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange) %>%
      mutate(Fast = Fast - y2$Medium - y3$Slow)

    # keeps <- grep(',1]', rownames(traces), invert=FALSE)
    # y6 <- as_tibble(traces[keeps,]) %>%
    #   select('2.5%', '25%', '50%', '75%', '97.5%') %>%
    #   gather('2.5%', '25%', '50%', '75%', '97.5%', key='pc', value='TF') %>%
    #   mutate(x=rep(calibxrange, times=5),
    #          xdate=as.Date(x-1, origin='1992-12-31')) %>%
    #   filter(x %in% calibxrange)

    keeps <- grep(",4]", rownames(traces), invert = FALSE)
    y4 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "TP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange) %>%
      group_by(x) %>%
      mutate(
        TPmed = median(TP),
        TPmin = min(TP),
        TPmax = max(TP),
        TPr = TP - TPmed
      )

    y4s <- y4 %>%
      select(x, xdate, pc, TP) %>%
      pivot_wider(names_from = pc, values_from = TP) %>%
      mutate(
        TPmed = `50%`,
        TPmin = `2.5%`,
        TPmax = `97.5%`,
        TF = y0$TF[y0$x == x]
      )

    adata$TP <- adata[[arun$chem1i]]/1000
    adata <- adata %>%
      group_by(x) %>%
      mutate(
        TPmed = median(y4$TP[y4$x == x]),
        TPr = TP - TPmed
      )

    keeps <- grep(",5]", rownames(traces), invert = FALSE)
    y5 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "TN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange) %>%
      group_by(x) %>%
      mutate(
        TNmed = median(TN),
        TNmin = min(TN),
        TNmax = max(TN),
        TNr = TN - TNmed
      )

    y5s <- y5 %>%
      select(x, xdate, pc, TN) %>%
      pivot_wider(names_from = pc, values_from = TN) %>%
      mutate(
        TNmed = `50%`,
        TNmin = `2.5%`,
        TNmax = `97.5%`,
        TF = y0$TF[y0$x == x]
      )

    adata$TN <- adata[[arun$chem2i]]/1000
    adata <- adata %>%
      group_by(x) %>%
      mutate(
        TNmed = median(y5$TN[y5$x == x]),
        TNr = TN - TNmed
      )

    keeps <- grep(",6]", rownames(traces), invert = FALSE)
    y6 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "fastTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",12]", rownames(traces), invert = FALSE)
    y12 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dfastTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    keeps <- grep(",7]", rownames(traces), invert = FALSE)
    y7 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "medTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",13]", rownames(traces), invert = FALSE)
    y13 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dmedTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    keeps <- grep(",8]", rownames(traces), invert = FALSE)
    y8 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "slowTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",14]", rownames(traces), invert = FALSE)
    y14 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dslowTP") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    keeps <- grep(",9]", rownames(traces), invert = FALSE)
    y9 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "fastTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",15]", rownames(traces), invert = FALSE)
    y15 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dfastTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    keeps <- grep(",10]", rownames(traces), invert = FALSE)
    y10 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "medTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",16]", rownames(traces), invert = FALSE)
    y16 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dmedTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    keeps <- grep(",11]", rownames(traces), invert = FALSE)
    y11 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "slowTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    keeps <- grep(",17]", rownames(traces), invert = FALSE)
    y17 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dslowTN") %>%
      mutate(
        x = rep(calibxrange, times = 5),
        xdate = as.Date(x-1, origin = refdate),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% calibxrange)
    
    ##### plot chem trace ####
    tci <- qt(.975, 7) # how many s.e. for 95% c.i.

    p4 <- ggplot() +
      labs(title = "", y = "DRP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      geom_ribbon(data = y4s, mapping = aes(x = xdate, ymin = TPmin - 0.02 * tci, ymax = TPmax + 0.02 * tci), colour = "lightgrey", fill = "lightgrey") +
      # geom_ribbon(data=y4s, mapping=aes(x=xdate, ymin=TPmin, ymax=TPmax), colour=tpcol[1], fill=tpcol[1]) +
      # geom_ribbon(data=y4s, mapping=aes(x=xdate, ymin=`25%`, ymax=`75%`), colour=tpcol[2], fill=tpcol[2]) +
      # geom_line(data=y4s, mapping=aes(x=xdate, y=TPmed), colour=tpcol[3]) +
      geom_line(data = y4, mapping = aes(x = xdate, y = TP, colour = pc)) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TP))

    p4r <- ggplot() +
      labs(title = "", y = "DRP residual " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      geom_ribbon(data = y4s, mapping = aes(x = xdate, ymin = TPmin - TPmed - 0.02 * tci, ymax = TPmax - TPmed + 0.02 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y4s, mapping = aes(x = xdate, ymin = TPmin - TPmed, ymax = TPmax - TPmed), colour = tpcol[1], fill = tpcol[1]) +
      geom_ribbon(data = y4s, mapping = aes(x = xdate, ymin = `25%` - TPmed, ymax = `75%` - TPmed), colour = tpcol[2], fill = tpcol[2]) +
      geom_line(data = y4s, mapping = aes(x = xdate, y = 0), colour = tpcol[3]) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TPr))

    robust <- function(formula, data, weights = weight) {
      rlm(formula, data, weights = weight, method = "M")
    }

    # y4s <- y4s %>% 
    #   filter(xdate >= dmy("01012005"))
    
    p4r2 <- ggplot() +
      labs(title = "", y = "DRP " ~ (mg ~ L^{-1}), x = "Model Median DRP " ~ (mg ~ L^{-1}), colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      # scale_x_date(limits=daterange, date_labels='%Y', date_breaks='1 year') +
      # scale_y_continuous(expand=c(0, 0), breaks=TPbreaks) +
      # coord_cartesian(ylim=TPlimits) +
      geom_ribbon(data = y4s, mapping = aes(x = TPmed, ymin = TPmin - 0.02 * tci, ymax = TPmax + 0.02 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y4s, mapping = aes(x = TPmed, ymin = TPmin, ymax = TPmax), colour = tpcol[1], fill = tpcol[1]) +
      geom_ribbon(data = y4s, mapping = aes(x = TPmed, ymin = `25%`, ymax = `75%`), colour = tpcol[2], fill = tpcol[2]) +
      geom_line(data = y4s, mapping = aes(x = TPmed, y = TPmed), colour = tpcol[3]) +
      # geom_smooth(data=adata, mapping=aes(x=TPmed, y=TP), linetype=2, colour="black", method="loess") +
      geom_point(data = adata, mapping = aes(x = TPmed, y = TP))

    flow_trans <- function() {
      # trans_new("flow", function(x) y=x, function(y) x=y)
      trans_new("flow", function(x) {
        y <- -1 / x
      }, function(y) {
        x <- -1 / y
      })
    }

    flow_breaks <- function(x) {
      r <- range(x)
      r <- r[r > 0 & is.finite(r)]
      round(rev(pretty(rev(r)^-1)^-1), 2)
      # b <- floor(log10(min(x))):ceiling(log10(max(x)))
      # b <- tcrossprod(10^b, c(1,2,5))
      # b <- sort(b)
      # b[b>=min(x) & b<=max(x)]
    }
    
    p4r3 <- ggplot() +
      labs(title = "", y = "DRP " ~ (mg ~ L^{-1}), x = "Flow " ~ (m^3 ~ s^{-1}), colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      # scale_x_date(limits=daterange, date_labels='%Y', date_breaks='1 year') +
      # scale_x_log10() +
      scale_x_continuous(limits = range(adata$TF), breaks = flow_breaks(adata$TF)) +
      # scale_x_continuous(limits = range(adata$TF[!is.na(adata$TP)])) +
      coord_trans(x = "flow") +
      # scale_y_continuous(expand=c(0, 0), breaks=TPbreaks) +
      # coord_cartesian(ylim=TPlimits) +
      geom_ribbon(data = y4s, mapping = aes(x = TF, ymin = TPmin - 0.02 * tci, ymax = TPmax + 0.02 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y4s, mapping = aes(x = TF, ymin = TPmin, ymax = TPmax), colour = tpcol[1], fill = tpcol[1]) +
      geom_ribbon(data = y4s, mapping = aes(x = TF, ymin = `25%`, ymax = `75%`), colour = tpcol[2], fill = tpcol[2]) +
      geom_line(data = y4s, mapping = aes(x = TF, y = TPmed), colour = tpcol[3]) +
      # geom_smooth(data=adata, mapping=aes(x=TPmed, y=TP), linetype=2, colour="black", method="loess") +
      geom_point(data = adata, mapping = aes(x = TF, y = TP))

    p6 <- ggplot() +
      labs(title = "", y = "fDRP" ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem1fast), shape = 1) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      # geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(high,1)), colour = wrtdscol, label = "95%") +
      geom_line(data = y6, mapping = aes(x = xdate, y = fastTP, colour = pc)) 
    p12 <- ggplot() +
      labs(title = "", y = "dfDRP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTPbreaks) +
      coord_cartesian(ylim = dTPlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y12, mapping = aes(x = xdate, y = dfastTP / 15, colour = pc))
    
    p7 <- ggplot() +
      labs(title = "", y = "mDRP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem1med), shape = 1) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(med,1)), colour = wrtdscol, label = "50%") +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y7, mapping = aes(x = xdate, y = medTP, colour = pc)) 
    p13 <- ggplot() +
      labs(title = "", y = "dmDRP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTPbreaks) +
      coord_cartesian(ylim = dTPlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y13, mapping = aes(x = xdate, y = dmedTP / 15, colour = pc))
    
    p8 <- ggplot() +
      labs(title = "", y = "sDRP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem1slow), shape = 1) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(low,1)), colour = wrtdscol, label = "5%") +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y8, mapping = aes(x = xdate, y = slowTP, colour = pc)) 
    p14 <- ggplot() +
      labs(title = "", y = "dsDRP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTPbreaks) +
      coord_cartesian(ylim = dTPlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y14, mapping = aes(x = xdate, y = dslowTP / 15, colour = pc))
    
    temp <- filter(adata, TP > 0.4)
    if (nrow(temp) > 0) {
      p4 <- p4 + # add labelled outliers
        geom_point(data = temp, mapping = aes(x = xdate, y = 0.37)) +
        geom_text(data = temp, mapping = aes(x = xdate, y = 0.37, label = as.character(TP)), nudge_x = 200)
      # p4r2 <- p4r2 + # add labelled outliers
      #   geom_point(data=temp, mapping=aes(x=TPmed, y=0.39), shape = 2)
      # geom_text(data=temp, mapping=aes(x=TPmed, y=0.37, label=as.character(TP)), nudge_y=0.01)
    }

    # testing
    # print(p4)

    #### trend analysis ####
    res4 <- y4 %>%
      filter(pc == "50%") %>%
      rename(TPmodel = TP) %>%
      left_join(adata, by = c("x", "xdate")) %>%
      drop_na() %>%
      mutate(TPres = TP - TPmodel) %>%
      mutate(
        xdate2 = decimal_date(xdate),
        month = month(xdate)
      )
    ggplot() +
      geom_point(data = res4, mapping = aes(x = xdate, y = TPres))
    if (sum(!is.na(res4$TP)) > 35) {
      res4fit <- rkt(date = res4$xdate2, y = res4$TPres, block = res4$month, correct = TRUE, rep ="m")
      print(paste(TPname, "res trend = ", res4fit$B / median(res4$TP, na.rm = TRUE), "p = ", res4fit$sl))
    } else if (sum(!is.na(res4$TP)) > 0) {
      res4fit <- rkt(date = res4$xdate2, y = res4$TPres)
      print(paste(TPname, "res trend = ", res4fit$B / median(res4$TP, na.rm = TRUE), "p = ", res4fit$sl))
    } else {
      print(paste("Not enough data for", TPname, "trend"))
    }

    p5 <- ggplot() +
      labs(title = "", y = "NO3-N " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = TNmin - 0.2 * tci, ymax = TNmax + 0.2 * tci), colour = "lightgrey", fill = "lightgrey") +
      # geom_ribbon(data=y5s, mapping=aes(x=xdate, ymin=TNmin, ymax=TNmax), colour=tncol[1], fill=tncol[1]) +
      # geom_ribbon(data=y5s, mapping=aes(x=xdate, ymin=`25%`, ymax=`75%`), colour=tncol[2], fill=tncol[2]) +
      # geom_line(data=y5s, mapping=aes(x=xdate, y=TNmed), colour=tncol[3]) +
      geom_line(data = y5, mapping = aes(x = xdate, y = TN, colour = pc)) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TN))

    # y5s <- y5s %>% 
    #   filter(xdate >= dmy("01012005"))
    
    p5r <- ggplot() +
      labs(title = "", y = "NO3-N residual " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = TNmin - TNmed - 0.2 * tci, ymax = TNmax - TNmed + 0.2 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = TNmin - TNmed, ymax = TNmax - TNmed), colour = tncol[1], fill = tncol[1]) +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = `25%` - TNmed, ymax = `75%` - TNmed), colour = tncol[2], fill = tncol[2]) +
      geom_line(data = y5s, mapping = aes(x = xdate, y = 0), colour = tncol[3]) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TNr))

    p5r2 <- ggplot() +
      labs(title = "", y = "NO3-N " ~ (mg ~ L^{-1}), x = "Model Median NO3-N" ~ (mg ~ L^{-1}), colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      # scale_x_date(limits=daterange, date_labels='%Y', date_breaks='1 year') +
      # scale_y_continuous(expand=c(0, 0), breaks=TNbreaks) +
      # coord_cartesian(ylim=TNlimits) +
      geom_ribbon(data = y5s, mapping = aes(x = TNmed, ymin = TNmin - 0.2 * tci, ymax = TNmax + 0.2 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y5s, mapping = aes(x = TNmed, ymin = TNmin, ymax = TNmax), colour = tncol[1], fill = tncol[1]) +
      geom_ribbon(data = y5s, mapping = aes(x = TNmed, ymin = `25%`, ymax = `75%`), colour = tncol[2], fill = tncol[2]) +
      geom_line(data = y5s, mapping = aes(x = TNmed, y = TNmed), colour = tncol[3]) +
      # geom_smooth(data=adata, mapping=aes(x=TNmed, y=TN), linetype=2, colour="black", method="loess") +
      geom_point(data = adata, mapping = aes(x = TNmed, y = TN))

    p5r3 <- ggplot() +
      labs(title = "", y = "NO3-N " ~ (mg ~ L^{-1}), x = "Flow " ~ (m^3 ~ s^{-1}), colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      # scale_x_date(limits=daterange, date_labels='%Y', date_breaks='1 year') +
      # scale_x_log10() +
      scale_x_continuous(limits = range(adata$TF), breaks = flow_breaks(adata$TF)) +
      # scale_x_continuous(limits = range(adata$TF[!is.na(adata$TN)])) +
      coord_trans(x = "flow") +
      # scale_y_continuous(expand=c(0, 0), breaks=TNbreaks) +
      # coord_cartesian(ylim=TNlimits) +
      geom_ribbon(data = y5s, mapping = aes(x = TF, ymin = TNmin - 0.2 * tci, ymax = TNmax + 0.2 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y5s, mapping = aes(x = TF, ymin = TNmin, ymax = TNmax), colour = tncol[1], fill = tncol[1]) +
      geom_ribbon(data = y5s, mapping = aes(x = TF, ymin = `25%`, ymax = `75%`), colour = tncol[2], fill = tncol[2]) +
      geom_line(data = y5s, mapping = aes(x = TF, y = TNmed), colour = tncol[3]) +
      # geom_smooth(data=adata, mapping=aes(x=TF, y=TN), linetype=2, colour="black", method="loess") +
      geom_point(data = adata, mapping = aes(x = TF, y = TN))

    p9 <- ggplot() +
      labs(title = "", y = "fNO3-N" ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem2fast), shape = 1) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      # geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(high,1)), colour = wrtdscol, label = "95%") +
      geom_line(data = y9, mapping = aes(x = xdate, y = fastTN, colour = pc)) 
    p15 <- ggplot() +
      labs(title = "", y = "dfNO3-N", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y15, mapping = aes(x = xdate, y = dfastTN / 15, colour = pc)) 
    
    p10 <- ggplot() +
      labs(title = "", y = "mNO3-N" ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem2med), shape = 1) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(med,1)), colour = wrtdscol, label = "50%") +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y10, mapping = aes(x = xdate, y = medTN, colour = pc)) 
    p16 <- ggplot() +
      labs(title = "", y = "dmNO3-N", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y16, mapping = aes(x = xdate, y = dmedTN / 15, colour = pc))
    
    p11 <- ggplot() +
      labs(title = "", y = "sNO3-N" ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      # geom_point(data = old_box, mapping = aes(x = xdate, y = chem2slow), shape = 1) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(low,1)), colour = wrtdscol, label = "5%") +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y11, mapping = aes(x = xdate, y = slowTN, colour = pc)) 
    p17 <- ggplot() +
      labs(title = "", y = "dsNO3-N", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y17, mapping = aes(x = xdate, y = dslowTN / 15, colour = pc))
    
    res5 <- y5 %>%
      filter(pc == "50%") %>%
      rename(TNmodel = TN) %>%
      left_join(adata, by = c("x", "xdate")) %>%
      # drop_na() %>%
      mutate(TNres = TN - TNmodel) %>%
      mutate(
        xdate2 = decimal_date(xdate),
        month = month(xdate)
      )
    ggplot() +
      geom_point(data = res5, mapping = aes(x = xdate, y = TNres))
    if (sum(!is.na(res4$TN)) > 35) {
      res5fit <- rkt(date = res5$xdate2, y = res5$TNres, block = res5$month, correct = TRUE, rep ="m")
      print(paste(TNname, "res trend = ", res5fit$B / median(res5$TN, na.rm = TRUE), "p = ", res5fit$sl))
    } else if (sum(!is.na(res4$TN)) > 0) {
      res5fit <- rkt(date = res5$xdate2, y = res5$TNres)
      print(paste(TNname, "res trend = ", res5fit$B / median(res5$TN, na.rm = TRUE), "p = ", res5fit$sl))
    } else {
      print(paste("Not enough data for", TNname, "trend"))
    }

    # fit plot
    # plotchem <- plot_grid(p4, p5, p0, nrow=3, align="v")
    # print(plotchem)
    # file_name <- paste(out_path, setname, '_chemtrace.png', sep="")
    # save_plot(file_name, plotchem, base_height=7, base_width=8)

    #### plot flow trace ####
    p1 <- ggplot() +
      labs(title = "", y = "Fast " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = fcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = flowbreaks) +
      coord_cartesian(ylim = flowlimits) +
      geom_line(data = y1, mapping = aes(x = xdate, y = Fast, colour = pc))

    p2 <- ggplot() +
      labs(title = "", y = "Med. " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = mcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = flowbreaks) +
      coord_cartesian(ylim = flowlimits) +
      geom_line(data = y2, mapping = aes(x = xdate, y = Medium, colour = pc))

    p3 <- ggplot() +
      labs(title = "", y = "Slow " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = scol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = flowbreaks) +
      coord_cartesian(ylim = flowlimits) +
      geom_line(data = y3, mapping = aes(x = xdate, y = Slow, colour = pc))

    print("WARNING: Simon commented out lines related to separation of hydrograph plot that weren't working")
    # y6$Slow <- filter(y3, pc == "50%")$Slow 
    # y6$Medium <- filter(y2, pc == "50%")$Medium
    # aah <- aah +
    #   geom_line(data = y6, mapping = aes(x = xdate, y = Medium + Slow), colour = "deepskyblue") + # add bach results
    #   geom_line(data = y6, mapping = aes(x = xdate, y = Slow), colour = "orchid") + # add bach results
    #   annotate(geom = "label", x = min(y6$xdate) + 120, y = max(TFlimits) * 0.8, label = setname, size = 4) # rewrite label

    # plotflow <- plot_grid(p1, p2, p3, nrow=3, align="v")
    # print(plotflow)
    # file_name <- paste(out_path, setname, '_flowtrace.png', sep="")
    # save_plot(file_name, plotflow, base_height=7, base_width=8)

    # p1 <- ggplot() +
    #   labs(title='', y='Fast '~(m^3~s^{-1}), x='', colour='Percentile') +
    #   theme_cowplot() +
    #   theme(legend.position='none', plot.margin=unit(c(0, 0.3, 0, 0), 'cm'), plot.title=element_blank()) +
    #   panel_border(colour='black') +
    #   scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
    #   scale_x_date(limits=daterange, date_labels='%Y', date_breaks='1 year') +
    #   scale_y_continuous(expand = c(0, 0), breaks=seq(0,10,2)) +
    #   coord_cartesian(ylim=c(0, 10)) +
    #   geom_line(data=y1, mapping=aes(x=xdate, y=Fast, colour=pc))

    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    plotall <- plot_grid(title, p4, p5, p0, p1, p2, p3, ncol = 1, rel_heights = c(0.3, 1, 1, 1, 1, 1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_alltrace.png", sep = "")
    save_plot(file_name, plotall, base_height = 10, base_width = 8)

    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    plotall <- plot_grid(title, p4, p5, ncol = 1, rel_heights = c(0.3, 1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_fittrace.png", sep = "")
    save_plot(file_name, plotall, base_height = 5, base_width = 8)
    
    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    plotall <- plot_grid(title, p4r, p5r, ncol = 1, rel_heights = c(0.3, 1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_residtime.png", sep = "")
    save_plot(file_name, plotall, base_height = 5, base_width = 8)
    
    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    blank <- ggdraw()
    plotall <- plot_grid(title, blank, p4r2, p5r2, ncol = 2, rel_heights = c(0.1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_residmedian.png", sep = "")
    save_plot(file_name, plotall, base_height = 5, base_width = 8)

    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    blank <- ggdraw()
    plotall <- plot_grid(title, blank, p4r3, p5r3, ncol = 2, rel_heights = c(0.1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_residflow.png", sep = "")
    save_plot(file_name, plotall, base_height = 5, base_width = 8)

    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    plotall <- plot_grid(title, p6, p7, p8, p9, p10, p11, ncol = 1, rel_heights = c(0.3, 1, 1, 1, 1, 1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_conctrace.png", sep = "")
    save_plot(file_name, plotall, base_height = 10, base_width = 8)

    title <- ggdraw() + draw_label(catchname, fontface = "bold")
    plotall <- plot_grid(title, p12, p13, p14, p15, p16, p17, ncol = 1, rel_heights = c(0.3, 1, 1, 1, 1, 1, 1), align = "v")
    # print(plotall)
    file_name <- paste(out_path, setname, "_dconctrace.png", sep = "")
    save_plot(file_name, plotall, base_height = 10, base_width = 8)
    
    # some useful information
    print(paste("Mean flow calib =", meanTFcalib))
    print(paste("Min fast flow =", min(y1$Fast)))
    print(paste("Max fast flow =", max(y1$Fast)))
    print(paste("Min medium flow =", min(y2$Medium)))
    print(paste("Max medium flow =", max(y2$Medium)))
    print(paste("Min slow flow =", min(y3$Slow)))
    print(paste("Max slow flow =", max(y3$Slow)))
  } # end if file exists

  ah[[arun$setseq]] <- aah # store hydrograph
} # end for loop

#### plot all hydrographs ####
if (nruns == 24) {
  blank <- ggplot(data.frame()) + geom_point()
  plotah <- plot_grid(
    blank, blank, blank,
    ah[[1]], ah[[2]], ah[[3]],
    ah[[4]], ah[[5]], ah[[6]],
    ah[[7]], ah[[8]], ah[[9]],
    ah[[10]], ah[[11]], ah[[12]],
    ah[[13]], ah[[14]], ah[[15]],
    ah[[16]], ah[[17]], ah[[18]],
    ah[[19]], ah[[20]], ah[[21]],
    ah[[22]], ah[[23]], ah[[24]],
    rel_heights = c(0.2, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 3, align = "v"
  )
  # print(plotah)
  file_name <- paste(out_path, "all_hydrotrace.png", sep = "")
  save_plot(file_name, plotah, base_height = 12, base_width = 12)
} else if (nruns == 8) {
  blank <- ggplot(data.frame()) + geom_point()
  plotah <- plot_grid(
    blank,
    ah[[1]], ah[[2]], ah[[3]],
    ah[[4]], ah[[5]], ah[[6]],
    ah[[7]], ah[[8]],
    rel_heights = c(0.2, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, align = "v"
  )
  # print(plotah)
  file_name <- paste(out_path, "all_hydrotrace.png", sep = "")
  save_plot(file_name, plotah, base_height = 12, base_width = 12)
}
