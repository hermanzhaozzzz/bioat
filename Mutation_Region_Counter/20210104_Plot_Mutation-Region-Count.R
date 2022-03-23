library(ggplot2)
library(tidyverse)
setwd("/Users/mac/mac_data/3.projects/25.2020-12_02_CBE_DetectSeq_paper_revise/MutationRegionCount.table/")
df <- read_tsv("293T-bat_EMX1-All-PD_rep1_max_distance_count.tsv",  )

v_max_count <- c()

for (row in 1:nrow(df)) {
  df_line <- df[row,]
  print(row)
  if(df_line$count_C2T==0 & df_line$count_G2A==0) {
    v_max_count <- c(v_max_count, 0)
  } else if (df_line$count_C2T!=0 & df_line$count_G2A==0) {
    v_max_count <- c(v_max_count, df_line$count_C2T %>% sum())
  } else if (df_line$count_C2T==0 & df_line$count_G2A!=0) {
    v_max_count <- c(v_max_count, df_line$count_G2A %>% sum())
  } else {
    v_max_count <- c(v_max_count, max(df_line$count_C2T %>% sum(), df_line$count_G2A %>% sum()))
  }
}

df$max_count <- v_max_count




ggplot(mtcars,aes(x =mpg))+
  geom_histogram(aes(y=..density..), # 纵坐标是密度。类似也可以将纵坐标设置为频数(count)
                 color="#88ada6", fill="#fffbf0", # 边框与填充色，可以不设置
                 alpha=.25,  # 透明度，可以不设置
                 binwidth = 2, # 柱子的宽度。类似得也可以设置柱子的个数，如bins = 30
                 center = 0) # 柱子与对应横坐标的相对位置。0是指居中对齐。1是指对应数字在柱子的右侧边线。可以不设置




