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

##### read in data ####
source("read_data9.r")
nruns <- nrow(runlist)

# choose rows to run
rows <- 1:nruns # to run all sites
ah <- vector("list", nruns) # all hydrographs
# rows <- c(20,21) # to run a subset of sites, useful for testing
# i <- rows[[10]]

#### loop through data sets ####
i <- rows[1]
for (i in rows) {

  # set axis options for each set
  if (i %in% c(1, 9, 17)) { # Tahuna
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 20)
    flowbreaks <- seq(0, 20, 5)
  } else if (i %in% c(2, 10, 18)) { # Otama
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 5)
    flowbreaks <- seq(0, 5, 1)
  } else if (i %in% c(3, 11, 19)) { # Waiotapu
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 10)
    flowbreaks <- seq(0, 10, 2)
  } else if (i %in% c(4, 12, 20)) { # Pokai
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 20)
    flowbreaks <- seq(0, 20, 5)
  } else if (i %in% c(5, 13, 21)) { # Piako
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 8)
    flowbreaks <- seq(0, 8, 2)
  } else if (i %in% c(6, 14, 22)) { # Waitoa
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 8)
    flowbreaks <- seq(0, 8, 2)
  } else if (i %in% c(7, 15, 23)) { # Waihou
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 80)
    flowbreaks <- seq(0, 80, 20)
  } else if (i %in% c(8, 16, 24)) { # Puniu
    TPlimits <- c(0, 0.4)
    TPbreaks <- seq(0, 0.8, 0.2)
    TNlimits <- c(0, 6)
    TNbreaks <- seq(0, 6, 2)
    flowlimits <- c(0, 60)
    flowbreaks <- seq(0, 60, 20)
  } else {
    stop()
  }

  dTPlimits <- c(-0.1, 0.1)
  dTPbreaks <- seq(-0.1, 0.1, 0.05)
  dTNlimits <- c(-1, 1)
  dTNbreaks <- seq(-1, 1, 0.5)
  
  # rough estimate for total flow (TF)
  TFlimits <- flowlimits * 1
  TFbreaks <- flowbreaks * 1

  data_file_name <- paste(runlist$catchfile[i], "_data.dat", sep = "")
  opt_file_name <- runlist$optfile[i]
  print(paste(data_file_name, "+", opt_file_name))

  # assemble run options/data into vectors (?)
  arun <- runlist[i, ]
  adata <- data[grep(pattern = data_file_name, x = data$file), ]
  aarea <- tibble(area = adata$area[1])
  aoptions <- options[opt_file_name, ]
  aalloptions <- cbind(arun, aoptions, aarea) # combine into one data table
  startcalib <- aalloptions$startcalib
  endcalib <- aalloptions$endcalib
  catchname <- aalloptions$catchname %>% print()
  setname <- aalloptions$setname
  keeps <- c("startrun", "startcalib", "endcalib", "startvalid", "endvalid")
  intoptions <- as.vector(t(aalloptions[keeps]))
  nintoptions <- length(intoptions)
  keeps <- c("chem1ae", "chem1re", "chem2ae", "chem2re", "area")
  realoptions <- as.vector(t(aalloptions[keeps]))
  nrealoptions <- length(realoptions)
  date <- as.vector(adata$date)
  ndate <- length(date)
  flow <- as.vector(adata$flow)
  TP <- as.vector(adata$TP)
  TN <- as.vector(adata$TN)
  meanTPcalib <- mean(TP[startcalib:endcalib], na.rm = TRUE)
  meanTNcalib <- mean(TN[startcalib:endcalib], na.rm = TRUE)
  meanTFcalib <- mean(flow[startcalib:endcalib])
  TP[is.na(TP)] <- -1 # stan doesn't understand NA
  TN[is.na(TN)] <- -1 # stan doesn't understand NA

  # add flow info
  adata <- adata %>%
    mutate(
      x = rep(1:ndate, times = 1),
      xdate = as.Date(date, origin = "1899-12-30"),
      TF = flow
    )

  # define date range
  fromdate <- as.Date(intoptions[2], origin = "2000-12-31")
  todate <- as.Date(intoptions[3], origin = "2000-12-31")
  daterange <- c(fromdate, todate + months(6))

  #### do hydrograph ####

  # get hydrograph
  y0 <- adata %>%
    mutate(pc = "data") %>%
    select(pc, TF, x, xdate) %>%
    filter(x %in% c(startcalib:endcalib)) %>%
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
    scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
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

    # Murdoch and Shanley trend analysis
    if (FALSE) {
      print("Mudoch and Shanley trend analysis")
      tdata <- adata %>%
        select(shortname, TF, TP, TN, x, xdate) %>%
        filter(x %in% c(startcalib:endcalib)) %>%
        mutate(
          day = day(xdate),
          month = month(xdate),
          year = year(xdate)
        )
      qflows <- quantile(tdata$TF, c(0.05, 0.50, 0.95))
      tyear <- tibble(year = unique(tdata$year)) %>%
        mutate(
          shortname = tdata$shortname[1],
          TP05 = NA_real_,
          TP50 = NA_real_,
          TP95 = NA_real_,
          TN05 = NA_real_,
          TN50 = NA_real_,
          TN95 = NA_real_
        )
      j <- 1
      for (j in 1:nrow(tyear)) {
        tdatai <- tdata %>%
          filter(year == tyear$year[j]) %>%
          drop_na()
        # modelTP <- loess(TP ~ TF, data = tdatai, control = loess.control(surface = "interpolate"))
        # modelTN <- loess(TN ~ TF, data = tdatai, control = loess.control(surface = "interpolate"))
        if (nrow(tdatai) > 3) {
          modelTP <- lm(TP ~ poly(TF, 3), data = tdatai)
          modelTN <- lm(TN ~ poly(TF, 3), data = tdatai)
          tdatai$modelTP <- predict(modelTP, tdatai)
          tdatai$modelTN <- predict(modelTN, tdatai)
          ggplot(tdatai) +
            geom_point(mapping = aes(x = TF, y = TP)) +
            geom_line(mapping = aes(x = TF, y = modelTP), colour = "red")
          ggplot(tdatai) +
            geom_point(mapping = aes(x = TF, y = TN)) +
            geom_line(mapping = aes(x = TF, y = modelTN), colour = "blue")
          tyear$TP05[j] <- predict(modelTP, data.frame(TF = qflows[1]))
          tyear$TP50[j] <- predict(modelTP, data.frame(TF = qflows[2]))
          tyear$TP95[j] <- predict(modelTP, data.frame(TF = qflows[3]))
          tyear$TN05[j] <- predict(modelTN, data.frame(TF = qflows[1]))
          tyear$TN50[j] <- predict(modelTN, data.frame(TF = qflows[2]))
          tyear$TN95[j] <- predict(modelTN, data.frame(TF = qflows[3]))
        }
      }
      ggplot(tyear) +
        labs(title = paste(tyear$shortname[1], "TP")) +
        geom_line(mapping = aes(x = year, y = TP05), colour = "blue") +
        geom_line(mapping = aes(x = year, y = TP50), colour = "green") +
        geom_line(mapping = aes(x = year, y = TP95), colour = "red")
      ggplot(tyear) +
        labs(title = paste(tyear$shortname[1], "TN")) +
        geom_line(mapping = aes(x = year, y = TN05), colour = "blue") +
        geom_line(mapping = aes(x = year, y = TN50), colour = "green") +
        geom_line(mapping = aes(x = year, y = TN95), colour = "red")
    }

    # USGS trend analysis
    if (TRUE) {
      print("USGS trend analysis")
      tdata <- adata %>%
        select(shortname, TF, TP, TN, x, xdate) %>%
        filter(x %in% c(startcalib:endcalib)) %>%
        mutate(
          day = day(xdate),
          month = month(xdate),
          year = year(xdate)
        )
      # convert my data to EGRET format
      Solute <- "TP"
      Solute <- "TN"
      wrtds <- list()
      for (Solute in c("TP", "TN")) {
        Daily <- data.frame(Date = tdata$xdate, Flow = tdata$TF)
        write_csv(Daily, "temp_egret_daily.csv")
        Daily <- readUserDaily(getwd(), "temp_egret_daily.csv", qUnit = 2, verbose = FALSE)
        Sample <- data.frame(Date = tdata$xdate, Remark = "", Conc = tdata[[Solute]]) %>% drop_na()
        write_csv(Sample, "temp_egret_sample.csv")
        Sample <- readUserSample(getwd(), "temp_egret_sample.csv")
        INFO <- data.frame(
          param.units = "mg/L",
          shortName = aalloptions$catchname,
          paramShortName = Solute,
          drainSqKm = aalloptions$area,
          constitAbbrev = Solute,
          staAbbrev = tdata$shortname[1],
          paStart = 1,
          paLong = 12,
          queryTime = Sys.time(),
          stringsAsFactors = FALSE
        )
        write_csv(INFO, "temp_egret_info.csv")
        INFO <- readUserInfo(getwd(), "temp_egret_info.csv", interactive = FALSE)
        eList <- mergeReport(INFO, Daily, Sample, verbose = FALSE)
        # plots and analysis
        eList <- modelEstimation(eList, verbose = FALSE, minNumObs = 60, minNumUncen = 30)
        if (tdata$shortname[1] %in% c("Ot", "Wp") && Solute %in% "TP"){
          # handle missing TP data, see EGRET manual p49
          eList <- blankTime(eList, "2004-12-21", "2011-10-25") # deletes poor data period
        }
        if (FALSE){
          multiPlotDataOverview(eList)
          plotConcTimeDaily(eList)
          plotConcPred(eList)
          plotResidPred(eList)
          plotResidQ(eList)
          plotResidTime(eList)
          boxResidMonth(eList)
          plotConcHist(eList)
          plotFluxHist(eList)
        }
        q1 <- quantile(eList$Daily$Q, 0.05)
        q2 <- quantile(eList$Daily$Q, 0.5)
        q3 <- quantile(eList$Daily$Q, 0.95)
        centerDate <- "07-01"
        yearStart <- eList$INFO$DecLow
        yearEnd <- eList$INFO$DecHigh
        # capture output of plot
        # print(Solute)
        wrtdscol <- "darkgrey"
        wrtdswt <- 0.5
        wrtds[[Solute]] <- plotConcTimeSmooth(eList, q1, q2, q3, centerDate, yearStart, yearEnd, 
                                              printValues = TRUE, 
                                              minNumObs = 60, 
                                              minNumUncen = 30, 
                                              windowY = 5) %>% # key!
          setNames(c("year", "low", "med", "high")) %>% 
          mutate(date = ymd(paste0(floor(year), centerDate))) 
        if (FALSE){
          fluxBiasMulti(eList)
          clevel <- pretty(range(eList$Sample$ConcAve))
          qlevel <- 10^pretty(log10(range(eList$Daily$Q)))
          plotContours(eList, yearStart, yearEnd, qBottom = qlevel[1], qTop = qlevel[length(qlevel)], contourLevels = clevel)
          plotDiffContours(eList, year0 = yearStart, year1 = yearEnd, qBottom = qlevel[1], qTop = qlevel[length(qlevel)], maxDiff = diff(range(eList$Sample$ConcAve)))
        }
      }
    }

    # read traces
    traces <- read_rds(file_name)

    # read quartiles from previous run to compare
    # old_path <- "../bach_constant/run - eckhardt_priors_narrow/"
    old_path <- paste0(out_path, "old_quartiles/")
    old_quartiles <- paste0(old_path, aalloptions$setname, "_quartiles.rds")
    old_quartiles <- str_replace(old_quartiles, "\\(x\\)", c("(a)", "(b)", "(c)"))
    old_fits <- vector("list", 3)
    old_fits[[1]] <- readRDS(old_quartiles[1]) %>%
      mutate(xdate = dmy("1-July-2004")) %>%
      slice(4:8)
    if (file.exists(old_quartiles[2])) {
      old_fits[[2]] <- readRDS(old_quartiles[2]) %>%
        mutate(xdate = dmy("1-July-2009")) %>%
        slice(4:8)
    } else {
      old_fits[[2]] <- readRDS(old_quartiles[1]) %>%
        replace_with_na_all(~TRUE) %>%
        mutate(xdate = dmy("1-July-2009")) %>%
        slice(4:8)
    }
    old_fits[[3]] <- readRDS(old_quartiles[3]) %>%
      mutate(xdate = dmy("1-July-2014")) %>%
      slice(4:8)
    old_box <- bind_rows(old_fits) %>% mutate(pc = rep(c("2.5%", "25%", "50%", "75%", "97.5%"), 3))

    # this is probably a better approach, use a single data frame for all variables
    # temp <- traces %>%
    #   select('2.5%', '25%', '50%', '75%', '97.5%') %>%
    #   mutate(name=rownames(traces)) %>%
    #   gather('2.5%', '25%', '50%', '75%', '97.5%', key='pc', value='val') %>%
    #   mutate(var=rep(c('Fast', 'Medium', 'Slow', 'TP', 'TN'), times=ndate*5),
    #          date=rep(1:ndate, each=5, times=5))

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
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))

    keeps <- grep(",2]", rownames(traces), invert = FALSE)
    y2 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "Medium") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib)) %>%
      mutate(Medium = Medium - y3$Slow)

    keeps <- grep(",1]", rownames(traces), invert = FALSE)
    y1 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "Fast") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib)) %>%
      mutate(Fast = Fast - y2$Medium - y3$Slow)

    # keeps <- grep(',1]', rownames(traces), invert=FALSE)
    # y6 <- as_tibble(traces[keeps,]) %>%
    #   select('2.5%', '25%', '50%', '75%', '97.5%') %>%
    #   gather('2.5%', '25%', '50%', '75%', '97.5%', key='pc', value='TF') %>%
    #   mutate(x=rep(1:ndate, times=5),
    #          xdate=as.Date(x, origin='2000-12-31')) %>%
    #   filter(x %in% c(startcalib:endcalib))

    keeps <- grep(",4]", rownames(traces), invert = FALSE)
    y4 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "TP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib)) %>%
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
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib)) %>%
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
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",12]", rownames(traces), invert = FALSE)
    y12 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dfastTP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    keeps <- grep(",7]", rownames(traces), invert = FALSE)
    y7 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "medTP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",13]", rownames(traces), invert = FALSE)
    y13 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dmedTP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    keeps <- grep(",8]", rownames(traces), invert = FALSE)
    y8 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "slowTP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",14]", rownames(traces), invert = FALSE)
    y14 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dslowTP") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    keeps <- grep(",9]", rownames(traces), invert = FALSE)
    y9 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "fastTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",15]", rownames(traces), invert = FALSE)
    y15 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dfastTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    keeps <- grep(",10]", rownames(traces), invert = FALSE)
    y10 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "medTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",16]", rownames(traces), invert = FALSE)
    y16 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dmedTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    keeps <- grep(",11]", rownames(traces), invert = FALSE)
    y11 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "slowTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    keeps <- grep(",17]", rownames(traces), invert = FALSE)
    y17 <- as_tibble(traces[keeps, ]) %>%
      select("2.5%", "25%", "50%", "75%", "97.5%") %>%
      gather("2.5%", "25%", "50%", "75%", "97.5%", key = "pc", value = "dslowTN") %>%
      mutate(
        x = rep(1:ndate, times = 5),
        xdate = as.Date(x, origin = "2000-12-31"),
        pc = factor(pc, levels = c("2.5%", "25%", "50%", "75%", "97.5%")[sort], ordered = TRUE)
      ) %>%
      filter(x %in% c(startcalib:endcalib))
    
    ##### plot chem trace ####
    tci <- qt(.975, 7) # how many s.e. for 95% c.i.

    p4 <- ggplot() +
      labs(title = "", y = "TP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      geom_ribbon(data = y4s, mapping = aes(x = xdate, ymin = TPmin - 0.02 * tci, ymax = TPmax + 0.02 * tci), colour = "lightgrey", fill = "lightgrey") +
      # geom_ribbon(data=y4s, mapping=aes(x=xdate, ymin=TPmin, ymax=TPmax), colour=tpcol[1], fill=tpcol[1]) +
      # geom_ribbon(data=y4s, mapping=aes(x=xdate, ymin=`25%`, ymax=`75%`), colour=tpcol[2], fill=tpcol[2]) +
      # geom_line(data=y4s, mapping=aes(x=xdate, y=TPmed), colour=tpcol[3]) +
      geom_line(data = y4, mapping = aes(x = xdate, y = TP, colour = pc)) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TP))

    p4r <- ggplot() +
      labs(title = "", y = "TP residual " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
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
      labs(title = "", y = "TP " ~ (mg ~ L^{-1}), x = "Model Median TP " ~ (mg ~ L^{-1}), colour = "Percentile") +
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
      labs(title = "", y = "TP " ~ (mg ~ L^{-1}), x = "Flow " ~ (m^3 ~ s^{-1}), colour = "Percentile") +
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
      labs(title = "", y = "fTP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem1fast), shape = 1) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(high,1)), colour = wrtdscol, label = "95%") +
      geom_line(data = y6, mapping = aes(x = xdate, y = fastTP, colour = pc)) 
    p12 <- ggplot() +
      labs(title = "", y = "dfTP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTPbreaks) +
      coord_cartesian(ylim = dTPlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y12, mapping = aes(x = xdate, y = dfastTP / 15, colour = pc))
    
    p7 <- ggplot() +
      labs(title = "", y = "mTP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem1med), shape = 1) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(med,1)), colour = wrtdscol, label = "50%") +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y7, mapping = aes(x = xdate, y = medTP, colour = pc)) 
    p13 <- ggplot() +
      labs(title = "", y = "dmTP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTPbreaks) +
      coord_cartesian(ylim = dTPlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y13, mapping = aes(x = xdate, y = dmedTP / 15, colour = pc))
    
    p8 <- ggplot() +
      labs(title = "", y = "sTP " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TPbreaks) +
      coord_cartesian(ylim = TPlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem1slow), shape = 1) +
      geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      geom_text(data = wrtds[["TP"]], mapping = aes(x = tail(date,1)+years(1), y = tail(low,1)), colour = wrtdscol, label = "5%") +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TP"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y8, mapping = aes(x = xdate, y = slowTP, colour = pc)) 
    p14 <- ggplot() +
      labs(title = "", y = "dsTP ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      # scale_colour_manual(values=c('orange','firebrick','darkred','firebrick','orange')) +
      scale_colour_manual(values = tpcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
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
      res4fit <- rkt(date = res4$xdate2, y = res4$TPres, block = res4$month, correct = TRUE)
      print(paste("TPres trend = ", res4fit$B / median(res4$TP, na.rm = TRUE), "p = ", res4fit$sl))
    } else if (sum(!is.na(res4$TP)) > 0) {
      res4fit <- rkt(date = res4$xdate2, y = res4$TPres)
      print(paste("TPres trend = ", res4fit$B / median(res4$TP, na.rm = TRUE), "p = ", res4fit$sl))
    } else {
      print("Not enough data for TP trend")
    }

    p5 <- ggplot() +
      labs(title = "", y = "TN " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
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
      labs(title = "", y = "TN residual " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = TNmin - TNmed - 0.2 * tci, ymax = TNmax - TNmed + 0.2 * tci), colour = "lightgrey", fill = "lightgrey") +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = TNmin - TNmed, ymax = TNmax - TNmed), colour = tncol[1], fill = tncol[1]) +
      geom_ribbon(data = y5s, mapping = aes(x = xdate, ymin = `25%` - TNmed, ymax = `75%` - TNmed), colour = tncol[2], fill = tncol[2]) +
      geom_line(data = y5s, mapping = aes(x = xdate, y = 0), colour = tncol[3]) +
      geom_point(data = adata, mapping = aes(x = xdate, y = TNr))

    p5r2 <- ggplot() +
      labs(title = "", y = "TN " ~ (mg ~ L^{-1}), x = "Model Median TN " ~ (mg ~ L^{-1}), colour = "Percentile") +
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
      labs(title = "", y = "TN " ~ (mg ~ L^{-1}), x = "Flow " ~ (m^3 ~ s^{-1}), colour = "Percentile") +
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
      labs(title = "", y = "fTN " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem2fast), shape = 1) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(high,1)), colour = wrtdscol, label = "95%") +
      geom_line(data = y9, mapping = aes(x = xdate, y = fastTN, colour = pc)) 
    p15 <- ggplot() +
      labs(title = "", y = "dfTN ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y15, mapping = aes(x = xdate, y = dfastTN / 15, colour = pc)) 
    
    p10 <- ggplot() +
      labs(title = "", y = "mTN " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem2med), shape = 1) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(med,1)), colour = wrtdscol, label = "50%") +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y10, mapping = aes(x = xdate, y = medTN, colour = pc)) 
    p16 <- ggplot() +
      labs(title = "", y = "dmTN ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y16, mapping = aes(x = xdate, y = dmedTN / 15, colour = pc))
    
    p11 <- ggplot() +
      labs(title = "", y = "sTN " ~ (mg ~ L^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = TNbreaks) +
      coord_cartesian(ylim = TNlimits) +
      geom_point(data = old_box, mapping = aes(x = xdate, y = chem2slow), shape = 1) +
      geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = low), colour = wrtdscol, linetype = 2, size = wrtdswt) +
      geom_text(data = wrtds[["TN"]], mapping = aes(x = tail(date,1)+years(1), y = tail(low,1)), colour = wrtdscol, label = "5%") +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = med), colour = wrtdscol, linetype = 5, size = wrtdswt) +
      # geom_line(data = wrtds[["TN"]], mapping = aes(x = date, y = high), colour = wrtdscol, linetype = 6, size = wrtdswt) +
      geom_line(data = y11, mapping = aes(x = xdate, y = slowTN, colour = pc)) 
    p17 <- ggplot() +
      labs(title = "", y = "dsTN ", x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = tncol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = dTNbreaks) +
      coord_cartesian(ylim = dTNlimits) +
      geom_hline(yintercept = 0) +
      geom_line(data = y17, mapping = aes(x = xdate, y = dslowTN / 15, colour = pc))
    
    res5 <- y5 %>%
      filter(pc == "50%") %>%
      rename(TNmodel = TN) %>%
      left_join(adata, by = c("x", "xdate")) %>%
      drop_na() %>%
      mutate(TNres = TN - TNmodel) %>%
      mutate(
        xdate2 = decimal_date(xdate),
        month = month(xdate)
      )
    ggplot() +
      geom_point(data = res5, mapping = aes(x = xdate, y = TNres))
    if (sum(!is.na(res4$TN)) > 35) {
      res5fit <- rkt(date = res5$xdate2, y = res5$TNres, block = res5$month, correct = TRUE)
      print(paste("TNres trend = ", res5fit$B / median(res5$TN, na.rm = TRUE), "p = ", res5fit$sl))
    } else if (sum(!is.na(res4$TN)) > 0) {
      res5fit <- rkt(date = res5$xdate2, y = res5$TNres)
      print(paste("TNres trend = ", res5fit$B / median(res5$TN, na.rm = TRUE), "p = ", res5fit$sl))
    } else {
      print("Not enough data for TN trend")
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
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = flowbreaks) +
      coord_cartesian(ylim = flowlimits) +
      geom_line(data = y1, mapping = aes(x = xdate, y = Fast, colour = pc))

    p2 <- ggplot() +
      labs(title = "", y = "Med. " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = mcol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
      scale_y_continuous(expand = c(0, 0), breaks = flowbreaks) +
      coord_cartesian(ylim = flowlimits) +
      geom_line(data = y2, mapping = aes(x = xdate, y = Medium, colour = pc))

    p3 <- ggplot() +
      labs(title = "", y = "Slow " ~ (m^3 ~ s^{-1}), x = "", colour = "Percentile") +
      theme_cowplot() +
      theme(legend.position = "none", plot.margin = unit(c(0, 0.3, 0, 0), "cm"), plot.title = element_blank()) +
      panel_border(colour = "black") +
      scale_colour_manual(values = scol[sort]) +
      scale_x_date(limits = daterange, date_labels = "%Y", date_breaks = "1 year") +
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

  ah[[aalloptions$setseq]] <- aah # store hydrograph
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
