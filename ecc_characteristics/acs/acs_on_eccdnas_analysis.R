library(data.table)
library(ggplot2)
library(scales)

df_large <- data.frame(fread('large_eccdnas_w_acs'))

large_eccdnas_w_acs_count <- nrow(df_large)

df_permute_large <- data.frame(fread('permute_large_eccs.txt'))

df_permute_large$freq <- df_permute_large$V1

p <- ggplot(df_permute_large, aes(x=freq)) + 
  geom_histogram(color="black", fill="white",bins=100) + theme_classic() +
  geom_vline(xintercept = large_eccdnas_w_acs_count, size=1.5) + 
  xlim(large_eccdnas_w_acs_count*0.5, large_eccdnas_w_acs_count*1.1) +
  labs(x="Regions containing ACS per permutation",
       y="Frequency")

ggsave("acs_on_large_eccs_permutation.pdf", plot = p, width = 6, height = 4)

df_large_srcount_w_acs <- data.frame(fread('large_eccdnas_w_acs_counted'))

df_large_srcount_w_acs$V5 <- rep("w_acs", nrow(df_large_srcount_w_acs))

colnames(df_large_srcount_w_acs) <- c('chrom', 'start','end','sr_count', 'acs')

df_large_srcount_w_acs$sr_count <- log10(df_large_srcount_w_acs$sr_count)
  
df_large_srcount_no_acs <- data.frame(fread('large_eccdnas_no_acs_counted'))

df_large_srcount_no_acs$V5 <- rep("no_acs", nrow(df_large_srcount_no_acs))

colnames(df_large_srcount_no_acs) <- c('chrom', 'start','end','sr_count', 'acs')

df_large_srcount_no_acs$sr_count <- log10(df_large_srcount_no_acs$sr_count)

df_large_srcount <- rbind(df_large_srcount_w_acs, df_large_srcount_no_acs)

p <- ggplot(df_large_srcount, aes(x=acs, y=sr_count))+geom_boxplot()+
  labs(x="",
       y="Log10 splitread count") + scale_x_discrete(labels=c("No ACS", "With ACS"))+
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8) , legend.title=element_text(size=8), legend.position = 'bottom', legend.text=element_text(size=6))+
  theme_classic()

p

ggsave("acs_on_large_eccs_srcount.pdf", plot = p, width = 2, height = 2)





df_micro <- data.frame(fread('micro_dnas_w_acs'))

micro_eccdnas_w_acs_count <- nrow(df_micro)

df_permute_micro <- data.frame(fread('permute_micro_dnas.txt'))

df_permute_micro$freq <- df_permute_micro$V1

p <- ggplot(df_permute_micro, aes(x=freq)) + 
  geom_histogram(color="black", fill="white",bins=100) + theme_classic() +
  geom_vline(xintercept = micro_eccdnas_w_acs_count, size=1.5) + 
  xlim(micro_eccdnas_w_acs_count*0.5, micro_eccdnas_w_acs_count*2) +
  labs(x="Regions containing ACS per permutation",
       y="Frequency")
p

ggsave("acs_on_micro_dnas_permutation.pdf", plot = p, width = 6, height = 4)

df_micro_srcount_w_acs <- data.frame(fread('micro_dnas_w_acs_counted'))

df_micro_srcount_w_acs$V5 <- rep("w_acs", nrow(df_micro_srcount_w_acs))

colnames(df_micro_srcount_w_acs) <- c('chrom', 'start','end','sr_count', 'acs')

df_micro_srcount_w_acs$sr_count <- log10(df_micro_srcount_w_acs$sr_count)

df_micro_srcount_no_acs <- data.frame(fread('micro_dnas_no_acs_counted'))

df_micro_srcount_no_acs$V5 <- rep("no_acs", nrow(df_micro_srcount_no_acs))

colnames(df_micro_srcount_no_acs) <- c('chrom', 'start','end','sr_count', 'acs')

df_micro_srcount_no_acs$sr_count <- log10(df_micro_srcount_no_acs$sr_count)

df_micro_srcount <- rbind(df_micro_srcount_w_acs, df_micro_srcount_no_acs)

p <- ggplot(df_micro_srcount, aes(x=acs, y=sr_count))+geom_boxplot()+
  labs(x="",
       y="Log10 splitread count") + scale_x_discrete(labels=c("No ACS", "With ACS"))+
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8) , legend.title=element_text(size=8), legend.position = 'bottom', legend.text=element_text(size=6))+
  theme_classic()

p

ggsave("acs_on_micro_dnas_srcount.pdf", plot = p, width = 2, height = 2)
