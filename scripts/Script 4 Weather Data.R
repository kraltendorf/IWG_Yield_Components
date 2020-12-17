# Project: IWG_Yield_Components
# Script 1 - Temperature and Precipitation Data
# Author: Kayla R. Altendorf
# Date: 12/15/2020

# load required packages
library("gridExtra")
library("tidyr")
library("tibble")
library("ggplot2")
library("dplyr")
library("cowplot")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG Yield Components/Scripts for GitHub/"

#### Step 1: Create a Dataframe of Events ####
# STP
event <- c("50% Anthesis", "50% Anthesis", "Harvest", "Harvest")
year <- as.factor(c(2017, 2018, 2017, 2018))
date <- c("2017-06-26", "2018-06-21", "2017-08-15", "2018-08-02")
day_of_year <- c(177, 172, 227, 214)
loc <- rep("STP", 4)

stp_events <- data.frame(event, year, date, day_of_year, loc)

# TLI
event <- c("50% Anthesis", "50% Anthesis", "Harvest", "Harvest")
year <- as.factor(c(2017, 2018, 2017, 2018))
date <- c("2017-06-10", "2018-06-06", "2017-07-25", "2018-07-24")
day_of_year <- c(161, 157, 206, 205)
loc <- rep("TLI", 4)

tli_events <- data.frame(event, year, date, day_of_year, loc)

# combine events together
events <- rbind(stp_events, tli_events)


#### Step 2: Prepare Weather Data ####
# this is from both STP and TLI
weather_dat <- read.csv(paste(dir, "data/Weather Data Abridged.csv", sep = ""), na.strings = c(".", NA, "")) %>% 
  mutate(tavg = ((as.numeric(TMAX) + as.numeric(TMIN)) / 2)) %>%
  separate(DATE, into = c("month", "day", "year"), sep = "/") %>% 
  mutate(year = as.numeric(year), 
         month = as.numeric(month), 
         day = as.numeric(day)) %>%
  arrange(STATION, year, month, day) %>%
  group_by(STATION, year) %>%
  mutate(day_of_year = row_number()) %>%
  select(STATION, month, day, year, day_of_year, PRCP_ACC, tavg)

# calculate ten year averages
ten_year_average <- weather_dat %>% 
  filter(year <= 16) %>%
  group_by(STATION, day_of_year) %>%
  summarise(PRCP_ACC_AVG = mean(PRCP_ACC, na.rm = T), TAVG = mean(tavg)) %>%
  mutate(year = "10-Year Average", 
         loc = case_when(STATION == "USC00218450" ~ "STP", 
                             STATION == "USW00003919" ~ "TLI"))

ten_year_average <- ten_year_average[,-1] # remove station column

experimental_years <- weather_dat %>% 
  filter(STATION == "USC00218450") %>%
  filter(year > 16) %>%
  group_by(year, day_of_year) %>%
  summarise(PRCP_ACC_AVG = mean(PRCP_ACC, na.rm = T), TAVG = mean(tavg)) %>%
  mutate(loc = "STP", 
         year = case_when(year == 16 ~ "2016", 
                          year == 17 ~ "2017", 
                          year == 18 ~ "2018")) %>%
  select(day_of_year, PRCP_ACC_AVG, TAVG, year, loc)

all_dat <- rbind(ten_year_average, experimental_years)

# since we have precision weather data taken right at the field, we can use TLI data for 2017 and 2018
tli16 <- read.csv(paste(dir, "data/", "TLI Weather Data 2016.csv", sep = ""), na.strings = c(".", "--"))
tli17 <- read.csv(paste(dir, "data/", "TLI Weather Data 2017.csv", sep = ""), na.strings = c(".", "--"))
tli18 <- read.csv(paste(dir, "data/", "TLI Weather Data 2018.csv", sep = ""), na.strings = c(".", "--"))
tli_dat <- bind_rows(tli16, tli17, tli18)

# group by date and get min and max per date
tli_dat <- tli_dat %>% 
  separate(col = date_time, into = c("date", "time"), sep = " ") %>%
  group_by(date) %>%
  mutate(temp_f = as.numeric(temp_f)) %>%
  summarise(TMAX = max(temp_f, na.rm = T), TMIN = min(temp_f, na.rm = T), PRCP = sum(rain_in, na.rm = T)) %>%
  ungroup() %>%
  mutate(TAVG = ((TMIN + TMAX) / 2)) %>%
  separate(col = date, into = c("month", "day", "year"), sep = "/") %>%
  arrange(year, as.numeric(month), as.numeric(day)) %>%
  mutate(loc = "TLI") %>%
  group_by(year) %>%
  mutate(day_of_year = row_number(), 
         PRCP_ACC_AVG = cumsum(PRCP)) %>%
  filter(year %in% c(17, 18)) %>%
  mutate(year = paste("20", year, sep = "")) %>%
  select(day_of_year, PRCP_ACC_AVG, TAVG, year, loc)

# change units to c and cm
all_dat <- rbind(all_dat, tli_dat) %>%
  mutate(TAVG = ((TAVG - 32) /1.8), 
         PRCP_ACC_AVG = PRCP_ACC_AVG * 2.54)



#### Step 3: Make the Figures ####
# temperature
font_size <- 25
temp <- ggplot(data = all_dat, aes(x = day_of_year, y = TAVG, color = year)) + 
  geom_line(size = 1) +
  facet_wrap(~loc, ncol = 1) + 
  geom_vline(aes(xintercept = day_of_year, linetype = event, color = year), size = 1, data = events, show.legend = TRUE) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  scale_color_manual(values = c("#616365", "#f0cc65", "#0072B2")) + 
  scale_x_continuous(name = "Day of Year") + 
  scale_y_continuous(name = "Average Daily Temperature (°C)", expand = c(0, 0), limits = c(-20, 40)) + 
  theme_minimal() + 
  theme(legend.text = element_text(size = font_size), 
        legend.position = "none",
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_text(size = font_size), 
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size))
  
# precipitation
precip <- ggplot(data = all_dat, aes(x = day_of_year, y = PRCP_ACC_AVG, color = year)) + 
  geom_line(size = 1) +
  facet_wrap(~loc, ncol = 1) + 
  geom_vline(aes(xintercept = day_of_year, linetype = event, color = year), size = 1, data = events, show.legend = TRUE) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  scale_color_manual(values = c("#616365", "#f0cc65", "#0072B2")) + 
  scale_x_continuous(name = "Day of Year") + 
  scale_y_continuous(name = "Accumulated Precipitation (cm)", expand = c(0, 0)) + 
  theme_minimal() + 
  theme(legend.text = element_text(size = font_size), 
        legend.position = "none", 
        legend.title = element_text(size = font_size), 
        strip.text.x = element_text(size = font_size), 
        axis.text.x = element_text(size = font_size), 
        axis.text.y = element_text(size = font_size), 
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size))


# extract legend (adding it back in the code above)
leg <- get_legend(precip)
# use cowplot to put them together
plots <- plot_grid(temp, precip, ncol = 2, labels = c("A", "B"), label_size = 25)
plot_grid(plots, leg, ncol = 1, rel_heights = c(1, .1))


