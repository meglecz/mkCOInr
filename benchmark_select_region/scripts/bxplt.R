library(ggplot2)

setwd("~/mkCOInr/benchmark_select_region/positive_mito_bait_select_region")

df<-read.table(file = file.path('count_trimmed_bait.tsv'), sep = '\t', header = TRUE)

df$bait <- paste(df$bait_taxlevel, df$N_seq_per_taxon, sep = "_")
df$identity <- paste(df$identity, "%", sep = " ")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
df$bait <- firstup(df$bait)

df$bait <- factor(df$bait , levels=c("Order_5", "Order_1", "Class_5", "Class_1", "Phylum_5", "Phylum_1"))
# One box per %identity
#png(file = "Sensitivity_per_bait_file.png")

p <- ggplot(df, aes(x=bait, y=trimmed_prop, fill=identity)) + 
  geom_boxplot(fill = "lightgray", color = "black") +
  facet_wrap(~identity) +
  labs(y="Sensitivity", x="Bait_type")
  p = p + ggtitle("Sensitivity for different bait sets and %identity thresholds") # for the main title
  p = p + theme_bw()
  p = p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
  p = p + theme(axis.text.y = element_text(size = 10))
  p = p + theme(axis.title.x = element_blank())
  p = p + theme(legend.position = "None")
  p = p + theme(plot.title = element_text(size = 10, hjust = 0.5))
p

ggsave("Sensitivity_per_bait_file.png", width = 12, height = 10, units = "cm")
# Save the file.
dev.off()



