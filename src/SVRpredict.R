#!/usr/bin/env Rscript

#version 1.19
suppressMessages(library(e1071))
suppressMessages(library(docopt))


'SVRpredict computes SVR score using an SVR model for 

Usage:
SVRpredict.R -i <input_sequences> -m <model> -o <output>
SVRpredict.R -i <input_sequences> -m <model>
SVRpredict.R (-h | --help)

Options:
-i  Input sequence file
-m  SVR model file
-o  Output file [default: prediction.txt]
-h  Show this screen.
--help    Show this screen.

Example of input sequence file:
CCCTATCGTGTAAACTATA
TTTTAGTCTATATCTGGTG
CACATACCAAGATAGTACT
GGAGATCCAGACAGCCGAT
' -> doc

arguments <- docopt(doc)
  if(arguments$o == F){arguments$output <- "prediction.txt"}


model <- readRDS(arguments$model)
mydf <- as.data.frame(read.table(arguments$input_sequences, stringsAsFactors = F, header = F), stringsAsFactors = F)
colnames(mydf) <- "DPR_sequence"


#formatting input for computing SVR score
dfsplit <- as.data.frame(strsplit(mydf[,1], ""), stringsAsFactors = F)
  dfsplit <- as.data.frame(t(dfsplit), stringsAsFactors = F)
  rownames(dfsplit) <- paste("sequence", 1:nrow(dfsplit), sep = "")

#print examples of sequences used
  cat("\n", "Sequence sample:", "\n")
  head(dfsplit)
  cat("\n")

#Initialize with polyA, polyC, polyG, polyT.
  polyNs <- as.data.frame(rbind(rep("A", ncol(dfsplit)),rep("C", ncol(dfsplit)),rep("G", ncol(dfsplit)),rep("T", ncol(dfsplit))), stringsAsFactor = F) 
  mydf_fact <- as.data.frame(rbind(polyNs, dfsplit), stringsAsFactor = F) 

#cConvert to 4-level factors    
  for(i in 1:ncol(mydf_fact)){
    mydf_fact[,i] <- as.factor(mydf_fact[,i])
  }
  mydf_fact <- mydf_fact[-(1:4),]

#predict SVR score
  pred <- as.data.frame(predict(model, mydf_fact), stringsAsFactor = F)
  colnames(pred) <- "SVR_score"
  outputdf <- cbind(mydf,round(pred, 2))

#write output to file
  write.table(outputdf, file = arguments$output, sep = '\t',row.names = F, quote = F)






