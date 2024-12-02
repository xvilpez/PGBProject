library(ggplot2)
library(dplyr)


#jaccard_data <- read.table("/home/leyao/Downloads/Tissues/jaccard_obox1_mef/jaccard_results.txt", header = FALSE, sep = "\t")
#jaccard_data <- read.table("/home/leyao/Downloads/Tissues/jaccard_myod1_mef/jaccard_results.txt", header = FALSE, sep = "\t")
#jaccard_data <- read.table("/home/leyao/Downloads/Tissues/jaccard_SRR396786/jaccard_results.txt", header = FALSE, sep = "\t")
jaccard_data <- read.table("/home/leyao/Downloads/Tissues/jaccard_myod1_all/jaccard_results.txt", header = FALSE, sep = "\t")

colnames(jaccard_data) <- c("Tissue", "Jaccard")

print(jaccard_data)

jaccard_data <- jaccard_data %>%
  arrange(Jaccard)

ggplot(jaccard_data, aes(x = reorder(Tissue, Jaccard), y = Jaccard, fill = Tissue)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "x",
    x = "Tissue", y = "Jaccard Index"
  ) +
  scale_fill_brewer(palette = "Set3") +
  coord_cartesian(ylim = c(0, 0.15))

