#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
library(readr)
library(reshape2)
library(tidyr)

##import csv file, split Sample Name column, remove no template control
input_file <- "data.csv" ### put the experiment name here
raw_data <- read_csv(input_file, skip=47, col_types = cols_only(
  'Sample Name' = col_character(),
  'CT' = col_number(),
  'Ct Mean'= col_number(),
  'Ct SD' = col_number(),
  'Target Name' = col_character()
))
raw_data <- as.data.frame(raw_data)
raw_data <- separate(raw_data, 'Sample Name', c('Genotype', 'Method'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- separate(raw_data, 'Target Name', c('Target', 'Name'), sep = ' ', extra = 'merge', remove = FALSE)
raw_data <- subset(raw_data, Genotype != 'no')
raw_data <- subset(raw_data, Target != 'YFP')

##calculate average from "pure" samples
pure_data <- subset(raw_data, Method == 'pure')
raw_data <- subset(raw_data, Method != 'pure')
raw_data <- raw_data[grep("dig", raw_data$Method),]
pure_average <- aggregate(pure_data, by = list(pure_data$Genotype), FUN = function(x) mean(as.numeric(as.character(x))))
pure_average[,3] <- pure_average[,1]
pure_average <- subset(pure_average, select = -c(Group.1))


##calculate dCT and remaining linear fraction
raw_data$CT <- raw_data$CT - pure_average$CT[match(raw_data$Genotype, pure_average$Genotype)]
raw_data[,10] <- 2 ^ -raw_data[,7] ## fraction

##set better column names, then assign time points to each genotype/method combination
raw_data[,11] <- NA
colnames(raw_data) <- c('samplename', 'genotype', 'method', 'target name', 'target', 'name', 'dct', 'ctmean', 'ctsd', 'fraction', 'timepoint')
for (row in 1:nrow(raw_data)) {
  genotype <- raw_data[row, 'genotype']
  method <- raw_data[row, 'method']
  if (genotype=='KA' | genotype=='GA')
  {
  if (method == 'dig 1') {
    raw_data[row,11] <- '48h'}
  else if(method == 'dig 2') {
    raw_data[row,11] <- '72h'}
  else if(method == 'dig 3') {
    raw_data[row,11] <- '96h'}
  else {next}
  }
  else
  {
    if (method == 'dig 1') {
      raw_data[row,11] <- '24h'}
    else if(method == 'dig 2') {
      raw_data[row,11] <- '48h'}
    else if(method == 'dig 3') {
      raw_data[row,11] <- '72h'}
    else {next}
  }
}

aggregate_raw_data <- aggregate(raw_data,by=list(raw_data$samplename,raw_data$genotype,raw_data$timepoint),FUN=mean)


##plot
ggplot(raw_data, aes(x = timepoint, y = fraction)) +
  geom_boxplot() + geom_point()

p <- ggplot(raw_data, aes(x=timepoint, y=fraction)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.5, width = 0.2) + theme_classic() + ylab('Linear DNA fraction') + xlab('') +
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8) , legend.title=element_text(size=8), legend.position = 'bottom', legend.text=element_text(size=6))


p

ggsave("qpcr_linear_dna_depletion_boxplot.pdf", plot = p, width = 2, height = 2)


# p <- ggplot(aggregate_raw_data, aes(x=Group.3, y=fraction)) + geom_jitter(size=0.5, width = 0.2) + theme_classic() + ylab('Linear DNA fraction') + xlab('') +
#   theme(axis.text=element_text(size=6),axis.title=element_text(size=8) , legend.title=element_text(size=8), legend.position = 'bottom', legend.text=element_text(size=6))
# 
# p


