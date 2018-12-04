library(data.table)
library(dplyr)
library(tidyr)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()


rm(list = ls())

path = "C:/Users/LakhanP/Desktop/HPC_stats"
setwd(path)

outFile = "HPC_usage_stats.xlsx"


userToLab = fread("fhs_userLabs.txt", sep = "\t", stringsAsFactors = F, header = T)


##########################################################################

## FHS_NORMAL
qNormal = fread("fhs_normal_queue_time.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(Submit = sub(pattern = "T\\d+:\\d+:\\d+", replacement = "", x = Submit, perl = T)) %>%
  mutate(Submit = as.Date(Submit, format = "%Y-%m-%d")) %>%
  mutate(year = format(Submit, "%Y"), month = month.abb[as.numeric(format(Submit, "%m"))]) %>%
  mutate(queue_time = as.numeric(queue_time))



userWiseNormalQ = qNormal %>% group_by(year, month, User) %>%
  summarise(NORMAL_qtime_total = as.numeric(sprintf("%.0f", sum(queue_time))), 
            NORMAL_qtime_avg = as.numeric(sprintf("%.0f", mean(queue_time))), 
            n = n()) %>%
  ungroup()

## FHS_NORMAL queue time: sum
normalqSumStats = userWiseNormalQ %>% 
  dplyr::select(month, User, NORMAL_qtime_total) %>%
  spread(month, NORMAL_qtime_total) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_NORMAL", type = "Total queue time") %>%
  dplyr::select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)

## FHS_NORMAL queue time: average
normalqMeanStats = userWiseNormalQ %>% 
  dplyr::select(month, User, NORMAL_qtime_avg) %>%
  spread(month, NORMAL_qtime_avg) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_NORMAL", type = "Average queue time") %>%
  dplyr::select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)

## FHS_NORMAL queue time: count
normalqNStats = userWiseNormalQ %>% 
  dplyr::select(month, User, n) %>%
  spread(month, n) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_NORMAL", type = "#in queue") %>%
  dplyr::select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)


FHS_NORMAL_qstat = bind_rows(normalqSumStats, normalqMeanStats, normalqNStats)


## FHS_NORMAL queue time: sum, average, count
normalqAllStats = userWiseNormalQ %>%
  gather(key = "variable", value = "value", c(NORMAL_qtime_total, NORMAL_qtime_avg, n)) %>%
  select(month, User, variable, value) %>%
  unite(temp, month, variable) %>%
  spread(temp, value) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_NORMAL") %>%
  select(Queue, Lab, User, Name, everything()) %>%
  arrange(Lab)



fwrite(x = FHS_NORMAL_qstat, file = "FHS_NORMAL_qstat.tab", sep = "\t", row.names = F, col.names = T, na = 0)
fwrite(x = normalqAllStats, file = "normalqAllStats.tab", sep = "\t", row.names = F, col.names = T, na = 0)



##########################################################################

## core usage stats
coresLong = fread("fhs_long.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  rename(FHS_LONG = `count(Partition)`)

coresNormal = fread("fhs_normal.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  rename(FHS_NORMAL = `count(Partition)`)


## per user stats
coreStatsByUser = coresLong %>% select(User, FHS_LONG) %>%
  full_join(y = coresNormal, by = ("User" = "User")) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  arrange(Lab) %>%
  mutate(FHS_LONG = ifelse(is.na(FHS_LONG),0,FHS_LONG), FHS_NORMAL = ifelse(is.na(FHS_NORMAL),0,FHS_NORMAL)) %>%
  select(User, Name, Lab, everything())

## per lab stats
coreStatsByLab = coreStatsByUser %>%
  group_by(Lab) %>%
  summarise(Users = n(), FHS_NORMAL = sum(FHS_NORMAL), FHS_LONG = sum(FHS_LONG)) %>%
  select(Lab, Users, everything())



fwrite(x = coreStatsByUser, file = "coreUsageStats_user.tab", sep = "\t", row.names = F, col.names = T)
fwrite(x = coreStatsByLab, file = "coreUsageStats_lab.tab", sep = "\t", row.names = F, col.names = T)




##########################################################################


## FHS_LONG queue time
qLong = fread("fhs_long_queue_time.csv", sep = ",", stringsAsFactors = F, header = T) %>%
  mutate(Submit = sub(pattern = "T\\d+:\\d+:\\d+", replacement = "", x = Submit, perl = T)) %>%
  mutate(Submit = as.Date(Submit, format = "%Y-%m-%d")) %>%
  mutate(year = format(Submit, "%Y"), month = month.abb[as.numeric(format(Submit, "%m"))]) %>%
  mutate(queue_time = as.numeric(queue_time))


userWiseLongQ = qLong %>% group_by(year, month, User) %>%
  summarise(LONG_qtime_total = as.numeric(sprintf("%.0f", sum(queue_time))), 
            LONG_qtime_avg = as.numeric(sprintf("%.0f", mean(queue_time))), 
            n = n()) %>%
  ungroup()


## FHS_LONG queue time: sum
longqSumStats = userWiseLongQ %>% select(month, User, LONG_qtime_total) %>%
  spread(month, LONG_qtime_total) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_LONG", type = "Total queue time") %>%
  select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)

## FHS_LONG queue time: average
longqMeanStats = userWiseLongQ %>% select(month, User, LONG_qtime_avg) %>%
  spread(month, LONG_qtime_avg) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_LONG", type = "Average queue time") %>%
  select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)

## FHS_LONG queue time: count
longqNStats = userWiseLongQ %>% select(month, User, n) %>%
  spread(month, n) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_LONG", type = "#in queue") %>%
  select(Queue, type, Lab, User, Name, everything()) %>%
  arrange(Lab)


FHS_LONG_qstat = bind_rows(longqSumStats, longqMeanStats, longqNStats)


## FHS_LONG queue time: sum, average, count
longqAllStats = userWiseLongQ %>%
  gather(key = "variable", value = "value", c(LONG_qtime_total, LONG_qtime_avg, n)) %>%
  select(month, User, variable, value) %>%
  unite(temp, month, variable) %>%
  spread(temp, value) %>%
  left_join(y = userToLab, by = c("User" = "ID")) %>%
  mutate(Queue = "FHS_LONG") %>%
  select(Queue, Lab, User, Name, everything()) %>%
  arrange(Lab)


fwrite(x = FHS_LONG_qstat, file = "FHS_LONG_qstat.tab", sep = "\t", row.names = F, col.names = T, na = 0)
fwrite(x = longqAllStats, file = "longqAllStats.tab", sep = "\t", row.names = F, col.names = T, na = 0)

