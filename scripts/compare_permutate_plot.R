#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
datFile <- args[1]
figureTitle <- paste('Folding energy for', args[2])

## data preprocess
dat <- read.csv(datFile, header=FALSE, sep='\t', stringsAsFactors = FALSE)
colnames(dat) <- c('Sequence', 'Energy', 'Length')
dat$seqNum <- sapply(dat$Sequence, function(x){
    names <- unlist(strsplit(x, '_'))
    if (length(names)==1){
        return(0)
    }else{
        return(as.numeric(names[-1]))
    }
})

### combine random sequences
dat$Type <- sapply(dat$Sequence, function(x){
    if (x=='Wildtype'){
        return('WildType')
    }else{
        return('Random')
    }     
})

### get minimum energy for each permutes
maxEdat <- dat %>%
    group_by(Sequence) %>%
    summarise(maxE=min(Energy))
maxEdat <- as.data.frame(maxEdat)
rownames(maxEdat) <- maxEdat$Sequence
dat$maxEnergy <- sapply(dat$Sequence, function(x){
    return(maxEdat[x,'maxE'])
})

## plot: separate sequences
ggplot(dat, aes(x=reorder(Sequence, seqNum), y=Energy, fill=maxEnergy)) + 
    geom_boxplot() +
    theme_bw() + 
    ggtitle(figureTitle) +
    xlab('Sequence ID') +
    ylab('Energy (kcal/mol)') +
    scale_fill_gradient(low='blue', high='red', name='MFE(kcal/mol)') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste(args[2], 'energy_separate.pdf', sep='_'), height = 4, width = 6, dpi = 300)

## plot: combine random
ggplot(dat, aes(x=Type, y=Energy, fill=Type)) +
    geom_violin() +
    geom_boxplot(width=0.1, fill='white') + 
    theme_bw() + 
    theme(legend.position = "none") +
    ggtitle(figureTitle) +
    ylab('Energy (kcal/mol)') +
    xlab('Sequence') +
    coord_flip()
ggsave(paste(args[2], 'energy_combined.pdf', sep='_'), height = 4, width = 5, dpi = 300)

