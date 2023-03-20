#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
NarnaFile <- args[1]
MitoFile <- args[2]

NarnaDat <- read.csv(NarnaFile, sep='\t', header=FALSE, stringsAsFactors = FALSE)
MitoDat <- read.csv(NarnaFile,  sep='\t', header=FALSE, stringsAsFactors = FALSE)

cols <- c('family', 'ID', 
         'seqLen', 'gcContent', 
         'wtMean', 'wtSd', 'wtNumStructure', 'wtMinEnergy', 
         'rdMean', 'rdSd', 'rdNumStructure', 'rdMinEnergy' )
colnames(NarnaDat) <- cols
colnames(MitoDat) <- cols

dat <- rbind(NarnaDat, MitoDat)

ggplot(dat) +
    geom_point(aes(x=wtMean, y=rdMean, color=gcContent, shape=family), size=2) +
    geom_abline(slope = 1, color="grey", linetype="dashed") +
    scale_color_gradient(low="blue", high="red", name='GC content') +
    scale_shape_discrete(labels=c('Mitovirus', 'Narnavirus'), name='Family') +
    ggtitle('Folding Energy to GC Content') +
    xlab('Folding Energy of Wildtype (kcal/mol)') +
    ylab('Folding Energy of Permutation (kcal/mol)') +
    xlim(-1550, -250) +
    ylim(-1550, -250) +
    theme_bw()

ggsave(filename = 'energy2gc.pdf', dpi = 300, width = 5.5, height = 4)