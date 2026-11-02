library(tidyverse)
library(ggridges)

files <- list.files(path = ".",
           pattern = "\\.xl\\.CPMnorm\\.bedgraph$",
           all.files = FALSE,
           full.names = TRUE,
           recursive = FALSE)

for (i in 1:length(files)){
  if (i == 1){
  df <- read.delim(files[i],
                       col.names = c("chr", "xlink_1", "xlink2", "count")) %>%
      mutate(count = abs(count),
             sample = str_remove(basename(files[i]), ".ncrna.xl.CPMnorm.bedgraph"))

  }

  else{
    df_2 <- read.delim(files[i],
                       col.names = c("chr", "xlink_1", "xlink2", "count")) %>%
      mutate(count = abs(count),
             sample = str_remove(basename(files[i]), ".ncrna.xl.CPMnorm.bedgraph"))

    df <- rbind(df, df_2)
  }
}


plot <- ggplot(df, aes(x = xlink_1, y = sample, weight = count)) +
  annotate("rect", xmin = 0,     xmax = 3656,  ymin = -Inf, ymax = Inf, fill = "#DEDDDD") +
  annotate("rect", xmin = 3657,     xmax = 5527,  ymin = -Inf, ymax = Inf, fill = "#E4F3F5") +
  annotate("rect", xmin = 5528,     xmax = 6622,  ymin = -Inf, ymax = Inf, fill = "#DEDDDD") +
  annotate("rect", xmin = 6623,     xmax = 6779,  ymin = -Inf, ymax = Inf, fill = "#E4F3F5") +
  annotate("rect", xmin = 6780,     xmax = 7934,  ymin = -Inf, ymax = Inf, fill = "#DEDDDD") +
  annotate("rect", xmin = 7935,     xmax = 12969,  ymin = -Inf, ymax = Inf, fill = "#E4F3F5") +
  annotate("rect", xmin = 12970,     xmax = 13314,  ymin = -Inf, ymax = Inf, fill = "#DEDDDD") +
  annotate("rect", xmin = 13315,     xmax = 42999,  ymin = -Inf, ymax = Inf, fill = "#DEDDDD") +
  geom_density_ridges(scale = 1, alpha = 0.8, bandwidth = 20) +
  theme_minimal() +
  ylab("Sample") +
  xlab("Crosslink position")


ggsave("rDNA_summary.png",
       plot = plot,
       dpi = 450,
       units = "cm")
