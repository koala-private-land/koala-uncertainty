library(tidyverse)
library(sf)
library(future.apply)
library(stringr)

# Load habitat suitability indices 

climateProjList = c("Avg", "CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3")
timestepList = c("t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7")
fmt <- "kitl_prop_%s_%s"

filenames <- lapply(climateProjList, function(cp) lapply(timestepList, function(t) sprintf(fmt, cp, t))) %>%
  unlist()

plan(multisession)
dbf <- future_lapply(filenames, function(f) {
  print(f)
  read_sf(dsn = "data/kitl_prop_unmasked.gdb", layer = f)
})

names(dbf) <- substr(filenames, 11, 25)

dbf_bind <- dbf %>%
  bind_rows(.id = "file") %>%
  mutate(timestep = str_sub(file, -2)) %>%
  mutate(climate_model = str_sub(file, 1, str_length(file)-3))

dbf_wide <- dbf_bind %>%
  select(-file, -AREA, -COUNT) %>%
  pivot_wider(names_from = timestep, values_from = MEAN)

# Koala habitat sizes in each property
khab_thres <- read_csv('data/khab_prop.csv') %>%
  rename(NewPropID = NEWPROPID) %>%
  mutate(khab_prop = SUM/COUNT)
