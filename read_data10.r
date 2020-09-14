# read data files into data frames
library(lubridate)
library(janitor)
library(readxl)
library(magrittr)
library(zoo)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

print("Reading Whatawhata data")

# avoid warning
as_numeric <- function(x, default = NA_real_) {
  suppressWarnings(dplyr::if_else(is.na(as.numeric(x)), default, as.numeric(x)))
}

# load run control file (put it in run dir so it gets saved with output)
file_name <- paste(out_path, "bachrunlist.tsv", sep = "")
runlist <- read_tsv(file = file_name, col_types = cols())

# add catchment names and ordering information
catch <- c("kiripaka" = 1, "mangaotama" = 2, "whakakai" = 3)
catchi <- catch[runlist$catchfile]
runlist <- runlist %>%
  mutate(
    catchname = c("Kiripaka", "Mangaotama", "Whakakai")[catchi],
    shortname = c("Ki", "Ma", "Wh")[catchi],
    catchseq = c(1, 2, 3)[catchi],
    setname = shortname,
    setseq = 1:nrow(runlist)
  )

# load conc data files
concdata <- vector("list", length(runlist$datafile))
i <- 1
for (i in match(unique(runlist$datafile), runlist$datafile)){ # don't need to repeat catchments
  suppressMessages({
    hdrs <- read_excel(paste0(out_path, runlist$datafile[i])) %>% clean_names() %>% names()
    raw <- read_excel(paste0(out_path, runlist$datafile[i]), col_names = FALSE, skip = 2) %>% 
      set_colnames(hdrs) %>% 
      select(-any_of(c("time", "coliforms_mpn_100ml", "vss"))) # drop problem columns
  })
  if (is.character(raw$date)){
    raw$date <- as.POSIXct(mdy(raw$date))
  }
  if (is.character(raw$date_time)){
    raw$date_time <- mdy_hm(raw$date_time)
  }
  temp <- raw %>% 
    mutate(
      time = hour(date_time) + minute(date_time)/60,
      p_h = as_numeric(p_h),
      catch = runlist$catchfile[i],
      catchname = runlist$catchname[i]
    )
  print(paste("Interpolating", sum(is.na(temp$time)), "missing sampling times")) 
  temp2 <- temp %>% 
    arrange(date) %>% 
    mutate(
      time = na.approx(time), # interpolate
      date_time = if_else(is.na(date_time), date+time*3600, date_time),
    )
  concdata[[i]] <- temp2
}
concdata <- bind_rows(concdata) 

# check missing
ggplot() +
  geom_line(data = concdata, mapping = aes(x = date, y = tp, colour = catch)) +
  geom_jitter(data = concdata %>% filter(is.na(drp)), mapping = aes(x = date, y = 1, colour = catch)) 
ggplot() +
  geom_line(data = concdata, mapping = aes(x = date, y = tn, colour = catch)) +
  geom_jitter(data = concdata %>% filter(is.na(no3_n)), mapping = aes(x = date, y = 1, colour = catch)) 

# date range
date_range <- range(c(ymd("19950101"), as.Date(concdata$date)))
date_range <- c(ymd("19930101"), ymd("20191231"))

# load flow data files
raw <- read_csv(paste0(out_path, "Flow observations_Whatatwhata.csv"), skip = 4, col_types = cols()) %>% clean_names()
raw <- raw %>% 
  mutate(
    date = mdy(date)
  ) %>% 
  pivot_longer(cols = all_of(runlist$catchfile), names_to = "catch", values_to = "flow") %>% 
  arrange(catch, date) %>% 
  filter(date >= date_range[1] & date <= date_range[2])
print(paste("Interpolating", sum(is.na(raw$flow)), "missing flow values")) 
flowdata <- raw %>% 
  group_by(catch) %>% 
  arrange(catch, date) %>% 
  mutate(
    lnflow = log(flow),
    lnflow = na.approx(lnflow, na.rm = FALSE), # warning could go <= 0
    flow = exp(lnflow)
  ) %>% 
  ungroup()
if (sum(is.na(flowdata$flow)) > 0){
  print(paste("Still", sum(is.na(flowdata$flow)), "missing flow values")) 
  print(paste("Using linear regression to interpolate missing Whakakai flow")) 
  temp <- flowdata %>% 
    arrange(catch, date) %>% 
    filter(date < ymd("20000101")) %>% # use this period for model
    pivot_wider(id_cols = "date", names_from = "catch", values_from = "lnflow")
  linmod <- lm(whakakai ~ kiripaka + mangaotama, data = temp)
  # summary(linmod)
  temp$whakakai_mod <- predict(linmod, newdata = temp)
  ggplot() +
    geom_line(data = temp, mapping = aes(x = date, y = whakakai), colour = "blue") +
    geom_line(data = temp, mapping = aes(x = date, y = whakakai_mod), colour = "red")
  i <- which(is.na(temp$whakakai))
  j <- which(flowdata$catch == "whakakai" & is.na(flowdata$lnflow))
  stopifnot(all(temp$date[i] == flowdata$date[j]))
  flowdata$lnflow[j] <- temp$whakakai_mod[i]
  flowdata$flow <- exp(flowdata$lnflow)
}

# check missing
ggplot() +
  geom_line(data = flowdata, mapping = aes(x = date, y = flow, colour = catch)) +
  scale_y_log10()





