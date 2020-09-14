# read data files into data frames
library(lubridate)
library(janitor)

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
    setseq = catchseq
  )

# load observed data files
data <- tibble()
raw <- read_csv(paste0(out_path, "Flow observations_Whatatwhata.csv"), skip = 4, col_types = cols()) %>% clean_names()
raw <- raw %>% 
  mutate(
    date = dmy(date)
  )
# check missing
ggplot() +
  




















dir_path <- out_path
files <- tibble(fname = list.files(path = dir_path, pattern = "\\_data.dat$"), path = dir_path)
for (i in 1:nrow(files)) {
  file_name <- paste(files$path[i], files$fname[i], sep = "")
  temp <- read_tsv(file_name,
    na = "NaN", skip = 0, progress = FALSE, n_max = Inf, col_types = cols(), guess_max = 5800,
    col_names = c(
      "date", "flow",
      "Ca", "Mg", "K", "Na", "Zn", "HCO3", "SO4", "Cl",
      "Si", "EC", "DRP", "NDP", "TP", "NNN", "NH4", "TN"
    )
  )
  keeps <- c("date", "flow", "TP", "TN") # columns to keep
  splitperiods <- FALSE
  temp <- temp[keeps] %>%
    mutate(
      file = files$fname[i], # add file name
      area = c(45, 103, 432, 519, 208, 802, 298, 121)[i],
      shortname = c("Ot", "Pi", "Po", "Pu", "Ta", "Wh", "Wp", "Wt")[i],
      catchseq = c(4, 8, 3, 6, 5, 1, 2, 7)[i],
      flow2 = flow / area * 86.4, # mm
      period = case_when(
        !splitperiods ~ "02-16",
        date >= 37257 & date <= 39082 ~ "02-06",
        date >= 39083 & date <= 40908 ~ "07-11",
        date >= 40909 & date <= 42735 ~ "12-16",
        TRUE ~ NA_character_
      ),
      period2 = case_when(
        !splitperiods ~ "x",
        date >= 37257 & date <= 39082 ~ "a",
        date >= 39083 & date <= 40908 ~ "b",
        date >= 40909 & date <= 42735 ~ "c",
        TRUE ~ NA_character_
      ),
      periodseq = case_when(
        !splitperiods ~ 1,
        date >= 37257 & date <= 39082 ~ 1,
        date >= 39083 & date <= 40908 ~ 2,
        date >= 40909 & date <= 42735 ~ 3,
        TRUE ~ NA_real_
      ),
      setname = paste(shortname, "(", period2, ")", sep = ""),
      rdate = as.Date(date, origin = "1899-12-30"),
      year = lubridate::year(rdate),
      setseq = periodseq + (catchseq - 1) * (nrow(runlist) / 8)
    )
  data <- rbind(data, temp)
}
