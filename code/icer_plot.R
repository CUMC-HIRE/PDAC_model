library(ggplot2)
library(tidyverse)
library(dplyr)


setwd("/Users/mai2125/cumc.columbia.edu/Lauren, Brianna N. - Pancreatic cancer/model/dump/icer")
df <- read_csv("panc_icers_bc_v2.csv")
df$arm <- recode(df$arm, 
                 ag = "G-nP",
                 folf = "FOLFIRINOX",
                 `natural history` = "Natural History")

ce_df <- df[!is.na(as.numeric(as.character(df$icers))),]

ggplot()+
  geom_line(data=ce_df, aes(x=cost/1000, y=QALYs))+
  geom_point(data=df,aes(x=cost/1000, y=QALYs, color=arm, shape=arm), size=4)+
  xlab("Cost (in thousands USD)")+
  labs(color="Strategy", shape="Strategy")+
  xlim(100, 250)+
  ylim(0, 3)

  
  
