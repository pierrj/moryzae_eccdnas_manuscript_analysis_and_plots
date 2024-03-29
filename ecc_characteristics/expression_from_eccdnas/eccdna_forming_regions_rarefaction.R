library(data.table)
library(ggplot2)
library(scales)

df_splitreads <- data.frame(fread('all.ecc_caller_out.splitreads.bed'))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()

df <- data.frame(fread('regions_found_rarefaction'))

df$means <- rowMeans(df[1:100])

df$sd <- apply(subset(df, select = 1:100), 1, sd)


df$steps <- seq(from = step_count, to = total_count, by = step_count)


p <- ggplot(df, aes(steps, means)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd))+
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x="Sampled splitreads",
       y="Unique splitreads")

p

ggsave("eccdna_forming_regions_rarefaction.pdf", plot = p, width = 6, height = 4)




df_splitreads <- data.frame(fread('all.ecc_caller_out.details.nolowq.txt'))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()

df <- data.frame(fread('details_regions_found_rarefaction'))

df$means <- rowMeans(df[1:100])

df$sd <- apply(subset(df, select = 1:100), 1, sd)


df$steps <- seq(from = step_count, to = total_count, by = step_count)


p <- ggplot(df, aes(steps, means)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd))+
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x="Sampled splitreads",
       y="Unique splitreads")

p

ggsave("eccdna_forming_regions_rarefaction_details.pdf", plot = p, width = 6, height = 4)



file <- 'ltr'

df_splitreads <- data.frame(fread(paste(file, '_eccdnas_all', sep='')))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()



df <- data.frame(fread(paste(file, '_regions_found_rarefaction', sep='')))

df$means <- rowMeans(df[1:100])

df$sd <- apply(subset(df, select = 1:100), 1, sd)


df$steps <- seq(from = step_count, to = total_count, by = step_count)


p <- ggplot(df, aes(steps, means)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd))+
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x="Sampled splitreads",
       y="Unique splitreads")

p

ggsave(paste("eccdna_forming_regions_rarefaction_", file, '.pdf', sep=''), plot = p, width = 6, height = 4)




file <- 'large'

df_splitreads <- data.frame(fread(paste(file, '_eccdnas_all', sep='')))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()



df <- data.frame(fread(paste(file, '_regions_found_rarefaction', sep='')))

df$means <- rowMeans(df[1:100])

df$sd <- apply(subset(df, select = 1:100), 1, sd)


df$steps <- seq(from = step_count, to = total_count, by = step_count)


p <- ggplot(df, aes(steps, means)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd))+
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x="Sampled splitreads",
       y="Unique splitreads")

p

ggsave(paste("eccdna_forming_regions_rarefaction_", file, '.pdf', sep=''), plot = p, width = 6, height = 4)




file <- 'micro'

df_splitreads <- data.frame(fread(paste(file, '_dnas_all', sep='')))

total_count <- nrow(df_splitreads)

step_count <- total_count/10

rm("df_splitreads")

gc()



df <- data.frame(fread(paste(file, '_regions_found_rarefaction', sep='')))

df$means <- rowMeans(df[1:100])

df$sd <- apply(subset(df, select = 1:100), 1, sd)


df$steps <- seq(from = step_count, to = total_count, by = step_count)


p <- ggplot(df, aes(steps, means)) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin=means-2*sd,ymax=means+2*sd))+
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x="Sampled splitreads",
       y="Unique splitreads")

p

ggsave(paste("eccdna_forming_regions_rarefaction_", file, '.pdf', sep=''), plot = p, width = 6, height = 4)