library(data.table)
library(ggplot2)
library(scales)

df_splitreads <- data.frame(fread('all.ecc_caller_out.regions.bed'))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()

df_observed <- data.frame(fread('genes_found_rarefaction'))

df_expected <- data.frame(fread('random_found_rarefaction'))


df_observed$means <- rowMeans(df_observed[1:100])

df_observed$sd <- apply(subset(df_observed, select = 1:100), 1, sd)

df_observed$steps <- seq(from = step_count, to = total_count, by = step_count)


df_expected$means <- rowMeans(df_expected[1:100])

df_expected$sd <- apply(subset(df_expected, select = 1:100), 1, sd)

df_expected$steps <- seq(from = step_count, to = total_count, by = step_count)

gene_count = 12115

p <- ggplot(df_observed, aes(steps, means)) + geom_point(aes(color="red")) + geom_line(aes(color="red")) + 
#  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd), width=30000)+
  geom_line(data=df_expected, aes(steps,means, color="blue")) +
  geom_point(data=df_expected, aes(steps,means, color="blue") ) +
#  geom_errorbar(data=df_expected, aes(ymin=means-2*sd, ymax=means+2*sd), width=30000) +
  scale_x_continuous(labels = comma, expand = expand_scale()) +
  scale_y_continuous(labels = comma, expand = expand_scale()) +
  labs(x="Sampled splitreads",
       y="Genes found on eccDNAs") + coord_cartesian(ylim = c(0, 12115)) +
  scale_color_identity(name = "",
                       breaks = c("red", "blue"),
                       labels = c("Observed", "Shuffled"),
                       guide = "legend")+theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.title=element_text(size=8))

p

ggsave("never_found_genes_rarefaction_analysis.pdf", plot = p, width = 3.25, height = 2)
