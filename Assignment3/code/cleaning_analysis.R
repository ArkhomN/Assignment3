library(dplyr)
library(ggplot2)
library(readr)

raw_data <- read_tsv("../raw_data/PlateletHW.tsv")

print(colnames(raw_data))

min_adp <- min(raw_data$ADP, na.rm = TRUE)

if (min_adp <= 0) {
  raw_data <- raw_data %>%
    mutate(ADP = ADP + abs(min_adp) + 1)
}

ggplot(raw_data, aes(x = ADP)) +
  geom_boxplot() +
  ggtitle("Boxplot of Transformed ADP Levels")

Q1 <- quantile(raw_data$ADP, 0.25)
Q3 <- quantile(raw_data$ADP, 0.75)
IQR <- Q3 - Q1
outlier_thresholds <- c(Q1 - 1.5 * IQR, Q3 + 1.5 * IQR)

clean_data <- raw_data %>%
  filter(ADP >= outlier_thresholds[1] & ADP <= outlier_thresholds[2])

write_csv(clean_data, "../clean_data/PlateletHW_cleaned.csv")

results <- list()
snp_vars <- c("rs4244285", "rs4986893", "rs662")

for (snp in snp_vars) {
  model <- lm(ADP ~ get(snp) + AGE + SEX, data = clean_data)
  results[[snp]] <- summary(model)
  print(paste("Results for SNP:", snp))
  print(summary(model))
}

summary_text <- "Association analysis between SNPs and ADP-induced platelet aggregation levels:\n"
for (snp in snp_vars) {
  snp_results <- results[[snp]]
  summary_text <- paste0(summary_text, 
                         "SNP: ", snp, "\n",
                         "  Estimate: ", snp_results$coefficients[2, 1], "\n",
                         "  Std. Error: ", snp_results$coefficients[2, 2], "\n",
                         "  t-value: ", snp_results$coefficients[2, 3], "\n",
                         "  p-value: ", snp_results$coefficients[2, 4], "\n\n")
}

write(summary_text, file = "../clean_data/association_summary.txt")
