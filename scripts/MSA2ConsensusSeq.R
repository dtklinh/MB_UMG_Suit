#!/usr/bin/Rscript

## calculate consensus seq from msa for a taxon
# 1st parameter: file path to msa of a taxon
# 2nd parameter: file to save consensus seq

library(Biostrings)
#library(msa)
library(dplyr)

## functions
msaConsensusSequence.matrix2 <- function(x, type=c("Biostrings", "upperlower"),
                                         thresh=c(80, 20), ignoreGaps=FALSE,
                                         ...)
{
  type <- match.arg(type)
  
  if (is.null(rownames(x)) ||
      any(!(rownames(x) %in% c(LETTERS, "-", "+", "."))))
    stop("consensus matrix 'x' is not in proper format")
  
  sel <- match(c("+", "."), rownames(x))
  sel <- sel[which(!is.na(sel))]
  if (length(sel) > 0)
    x <- x[-sel, ]
  
  if (!is.numeric(thresh) || length(thresh) !=2 ||
      thresh[1] > 100 || thresh[2] < 0 || thresh[1] < thresh[2])
    stop("'thresh' must be a decreasing sequence of two numbers ",
         "between 0 and 100")
  
  thresh <- thresh / 100
  
  if (type == "Biostrings")
  {
    cs <- colSums(x)
    
    if (any(is.na(cs)))
    {
      res <- rep.int("#", ncol(x))
      
      sel <- which(!is.na(cs))
      
      if (length(sel) > 0)
      {
        sstr <- consensusString(x[, sel, drop=FALSE], ...)
        res[sel] <- unlist(strsplit(sstr, ""))
      }
      
      out <- paste(res, collapse="")
    }
    else
      out <- consensusString(x, ...)
  }
  else
  {
    if (ignoreGaps)
    {
      sel <- match("-", rownames(x))
      
      if (!is.na(sel))
        sel <- (1:nrow(x))[-sel]
      else
        sel <- (1:nrow(x))
      
      perColFunc <- function(y)
      {
        if (any(is.na(y)))
          char <- "#"
        else
        {
          y <- y[sel] / sum(y[sel])
          
          maxw <-which.max(y[sel])
          maxi <- y[sel[maxw]]
          
          char <- rownames(x)[sel[maxw]]
          
          if (maxi < thresh[1])
          {
            if (maxi >= thresh[2])
              char <- tolower(char)
            else
              char <- "."
          }
        }
        
        char
      }
    }
    else
    {
      perColFunc <- function(y)
      {
        if (any(is.na(y)))
          char <- "#"
        else
        {
          y <- y / sum(y)
          
          maxw <- which.max(y)
          maxi <- y[maxw]
          
          char <- rownames(x)[maxw]
          
          if (maxi < thresh[1])
          {
            if (maxi >= thresh[2])
              char <- tolower(char)
            else
              char <- "."
          }
        }
        
        char
      }
    }
    
    out <- paste(apply(x, 2, perColFunc), collapse="")
  }
  
  out
}
##--------
msaConsensusSequence.MultipleAlignment2 <- function(x, type=c("Biostrings", "upperlower"),
                                                    thresh=c(80, 20), ignoreGaps=FALSE,
                                                    ...)
{
  m <- consensusMatrix(x)
  return(msaConsensusSequence.matrix2(m,type,thresh, ignoreGaps,...))
}
### end functions --------------------------------------

#setMethod("msaConsensusSequence2", signature("matrix"), definition=msaConsensusSequence.matrix2)


#setMethod("msaConsensusSequence2", signature("MultipleAlignment"),
# function(x, ...)
# {
#   mat <- consensusMatrix(x)
#   msaConsensusSequence2(mat, ...)
# })

args <- commandArgs(trailingOnly = TRUE)

taxon_msa <- Biostrings::readDNAMultipleAlignment(args[1])
TaxID <- names(taxon_msa@unmasked)[1]

taxon_consensus_seq <- msaConsensusSequence.MultipleAlignment2(taxon_msa, type="upperlower", thresh=c(40, 25))

consen_seq_trim <- gsub('[-?.]', '', taxon_consensus_seq) %>% toupper()
consen_seq_DNAStr <- DNAStringSet(consen_seq_trim)
names(consen_seq_DNAStr) <- toString(TaxID)
writeXStringSet(consen_seq_DNAStr, args[2])




