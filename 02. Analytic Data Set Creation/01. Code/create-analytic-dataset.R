# Program Namne: create-analytic-dataset.R
# Author: Jacob Englert
# Purpose: Compile meteorology, health, and time data for an analysis of
# Alzheimer's disease and heatwaves in California.


# Load Packages -----------------------------------------------------------
library(tidyverse)


# Load Data ---------------------------------------------------------------

# Location data (2015, US, Zip-code level)
l_data <- sf::st_read(here::here('01. Data','Location Data','ESRI ZIP 2015/ESRI15USZIP5_POLY_WGS84.shp'))
l_data_ca <- filter(l_data, STATE == 'CA') # 30,541 US ZIP codes --> 1,707 CA ZIP codes

# Meteorology data (2005-2015, CA, Zip-code level)
m_data <- read_rds(here::here('01. Data','Exposure Data','Daymet_CA.rda')) |>
  filter(year(DATE) %in% 2005:2015)
max(table(m_data$ZCTA, m_data$DATE)) # No duplicate entries

# ED visit data (2005-2015, CA)
h_data <- haven::read_sas(here::here('01. Data','Health Data','ca.sas7bdat')) |>
  filter(AD == 1) # subset to only AD patients
h_data[h_data == ""] <- NA

# Filter to only CA patient zip codes
# h_data_ca <- filter(h_data, patzip %in% l_data$ZIP) # 644,633 visits --> 636,660 visits
h_data_ca <- h_data

# Time data
t_data <- haven::read_sas(here::here('01. Data','Time Data','time_93_20_short.sas7bdat')) |>
  filter(year(DATE) %in% 2005:2015) |>
  select(DATE, HOLIDAY_FO)



# Manually create heatwave indicator and moving averages ------------------

expos_df <- m_data |>
  mutate(Year = year(DATE),
         Month = month(DATE),
         DoW = weekdays(DATE),
         DoY = yday(DATE)) |>
  arrange(ZCTA, DATE) |>
  mutate(extreme_heat = ifelse(AVG > quantile(AVG, 0.95, TRUE), 1, 0),
         AVG3DMA = zoo::rollmeanr(AVG, 3, fill = NA),
         DP3DMA = zoo::rollmeanr(DP, 3, fill = NA),
         .by = ZCTA) |>
  mutate(HW = zoo::rollsumr(extreme_heat, 2, fill = NA),
         .by = ZCTA) |>
  mutate(HW = ifelse(HW == 2, 1, 0))
# AVG3DMA: 3-day moving average for average daily temperature
# DP3DMA: 3-day moving average for daily dew point
# HW: heatwave indicator (two or more days above the 95th percentile of zip-code specific average temperature)

rm(m_data)

# Visualize exposure data
expos_df |>
  filter(ZCTA %in% sample(expos_df$ZCTA, 25)) |>
  ggplot(aes(x = DATE, y = AVG, color = factor(HW))) +
  geom_point() +
  facet_wrap(~ZCTA)


# Create Case-Crossover Dataset -------------------------------------------

# Clean up health data (create binary race/sex/ethnicity indicators)
h_data_ca <- h_data_ca |>
  mutate(WHITE = ifelse(RACE == 0, 1, 0),
         BLACK = ifelse(RACE == 1, 1, 0),
         AIAN = ifelse(RACE == 2, 1, 0),
         API = ifelse(RACE == 3, 1, 0),
         RACEOTH = ifelse(RACE %in% c(2,4), 1, 0),
         across(WHITE:RACEOTH, ~ ifelse(is.na(RACE), NA, .)),
         FEMALE = as.numeric(sex),
         HISPANIC = as.numeric(ETH)) |>
  select(ADMDATE = date, ZCTA = patzip, AGE = age, FEMALE,
         WHITE, BLACK, API, RACEOTH, HISPANIC, AD, ADRD, AD1, ADRD1,
         CKD:STROKE1)


cco <- h_data_ca |>
  filter(ZCTA %in% expos_df$ZCTA) |> # only include patients with exposure information
  mutate(ID = row_number(),
         Year = year(ADMDATE),
         Month = month(ADMDATE),
         DoW = weekdays(ADMDATE)) |>
  left_join(expos_df, by = c('Year','Month','DoW','ZCTA'), 
            relationship = 'many-to-many') |>
  mutate(Case = ifelse(DATE == ADMDATE, 1, 0)) |>
  select(-ADMDATE) |>
  left_join(t_data, by = 'DATE') |>
  select(ID, AGE:STROKE1, STATE, ZCTA, DATE, Case, Year:DoW, DoY, AVG:AT, AVG3DMA, DP3DMA, HW, HOLIDAY = HOLIDAY_FO)

# Export
write_csv(cco, here::here('02. Analytic Data Set Creation','02. Final Analytic Data Set','ad-ca-cco.csv'))
