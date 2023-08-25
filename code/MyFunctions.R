## my functions

## check if a case_id satisfies the paired PT-STN
Check_paired_PT_STN <- function(case_id, Tab){
  lst_sample <- Tab[Tab$case_id==case_id, "Sample"]
  return(("PT" %in% lst_sample) & ("STN" %in% lst_sample))
}
## remove technical replicate of sample
Remove_TechReplicateSample <- function(CountTable){
  ## keep sample with higher value
  Sample_bcr <- lapply(names(CountTable), function(x){return(substr(x,1,16))}) %>% unlist()
  df <- data.frame(Whole_SampleBCR=names(CountTable), SampleBCR=Sample_bcr)
  df2 <- aggregate(df$Whole_SampleBCR, by=list(df$SampleBCR), FUN=max)
  keep_sample_bcr <- df2[,2]
  return(CountTable[,c(keep_sample_bcr)])
}
##----------------------
fun <- function(test, x) {test[which.max(ave(test, test, FUN = seq_along) == x)]}
##---------------------------------
ExtractTaxonomy <- function(TaxoString, m_rank="k"){
  pat <- "__\\s*(.*?)\\s*;"
  if(m_rank %in% c("k", "p", "c", "o", "f", "g")){
    pat <- paste0(m_rank, pat)
  } else if(m_rank=="s"){
    return(strsplit(TaxoString, "s__")[[1]] %>% tail(n=1))
  } else{
    print("Wrong input!")
    return(NA)
  }
  res <- str_match(TaxoString, pat)
  return(res[2])
}
##-------------------
Extract_FastaID <- function(AllSeqNames, SpeciesName){
  is_this_species <- grepl(paste0(" ", SpeciesName, " "), AllSeqNames, fixed = T)
  Lst_id <- lapply(AllSeqNames[is_this_species], function(x){strsplit(x, " ")[[1]][1]}) %>% unlist()
  return(Lst_id)
}

##---------------------------------------
## Process and merge interval to infer 16S sequence
Infer_16S_Seq <- function(IR, min.gap = 300){
  ## IR contains all intervals from blastn for a particular taxon
  ## the first interval of IR supposed to be the best hits
  ##---------
  ## only keep the interval same strand with the first one
  IR <- IR[IR@strand==IR[1]@strand,]
  IR_disjoin <- reduce(IR,min.gapwidth=min.gap)
  return(IR_disjoin[which.max(IRanges::width(IR_disjoin))])
}
##---------
Validate_16S_Seq <- function(RefGenome, best_hit_16S, F.nuc_genome){
  ## validate if the inferred sequence is similar to 16S rRNA from F. Nucleatum
  #best_hit_16S <- Biostrings::extractAt(RefGenome[names(RefGenome)==(seqnames(best_hit) %>% toString())], best_hit@ranges)[[1]]
  consen_seq <- msa::msaConsensusSequence(msa(c(F.nuc_genome, best_hit_16S), method = "ClustalOmega"))
  lf <- letterFrequency(Biostrings::BString(consen_seq), letters = c("C", "G", "A", "T", "?"))
  return(lf["?"]/nchar(consen_seq))
}
##----------------------------
## Extract 16S sequence from DNAStringSet
Extract_16S <- function(StrSet, TaxID){
  idx <- which.min(abs(width(StrSet[StrSet@metadata$TaxID==TaxID])-1500))
  return(StrSet[StrSet@metadata$TaxID==TaxID][idx])
}
##-------------------------------------
Complete_Infer_16S <- function(Path2RefGenomeFaFile, F.nuc_16S_Seq){
  #---------- this function is to infer the 16S sequence of a taxon given its genome
  #---------- Path2RefGenomeFaFile: path to denome file, eg. 123.fa
  #---------- F.nuc_16S_Seq: DNAStringSet of F. nucleatum (or other well known species) 16S rRNA sequence.
  ## Check if NCBI balst path is in the PATH
  if(!(grepl("ncbi-blast", Sys.getenv("PATH"), fixed = T))){
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, "/sybig/home/lda/Software/ncbi-blast-2.14.0+/bin", sep = ":"))
  }
  ## read ref genome
  ref_genome <- readDNAStringSet(Path2RefGenomeFaFile)
  ## prepare database
  makeblastdb(Path2RefGenomeFaFile, dbtype = "nucl")
  ## blast and predict
  db <- blast(Path2RefGenomeFaFile)
  res <- predict(db, F.nuc_16S_Seq[1], BLAST_args = "-task blastn")
  ## pick up the most likely sequence
  target_Acc <- res$SubjectID[1]
  res <- res[res$SubjectID==target_Acc,]
  ## process to irange
  res <- tibble(res)
  res <- res %>% dplyr::select(SubjectID, S.start, S.end, E)
  res <- res %>% mutate(i_start = case_when(S.start<S.end ~S.start, S.end<S.start ~S.end))
  res <- res %>% mutate(i_end = case_when(S.start<S.end ~ S.end, S.end < S.start ~ S.start))
  res <- res %>% mutate(strand = case_when(S.start<S.end ~ "+", S.end<S.start ~ "-"))
  res <- res %>% select(SubjectID, i_start, i_end, strand)
  #BlastRes <- BlastRes[,c(2:5,1)]
  colnames(res) <- c("seqnames", "start", "end", "strand")
  res <- res[, c(2,3,1,4)]
  res_iranges <- as_granges(res)
  best_hit <- Infer_16S_Seq(res_iranges, 300)
  is_seq <- grepl((GenomeInfoDb::seqnames(best_hit) %>% toString()), names(ref_genome), fixed = T)
  best_hit_16S <- Biostrings::extractAt(ref_genome[is_seq], best_hit@ranges)[[1]]
  ## check strand
  if((best_hit@strand %>% toString())=="-"){
    best_hit_16S <- Biostrings::reverseComplement(best_hit_16S)
  }
  return(list(DNAStrSet=best_hit_16S, irange=best_hit))
}
####------------------------
Distance_1_to_Set <- function(seq, StrSet){
  lst_res <- c()
  N <- length(StrSet)
  
  doParallel::registerDoParallel(12)
  lst_res <- foreach(i = 1:N, .combine=c) %dopar% {
    m <- msaClustalOmega(c(seq, StrSet[i])) %>% DNAStringSet()
    dist <- DECIPHER::DistanceMatrix(m)
    #lst_res <- c(lst_res, dist[1,2])
    dist[1,2]
  }
  return(lst_res)
}
###----------------------------
perColFunc <- function(y)
{
  if (any(is.na(y)))
    char <- "#"
  else
  {
    y <- y / length(y)
    
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


