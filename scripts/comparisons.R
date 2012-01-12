
stringsAsFactors = FALSE

qw <- function(...) {
  as.character(sys.call()[-1])
}

ns5 <- read.csv(file = "data/ns5_nearest_peak_to_gene_TSS.csv")

astro <- read.csv(file = "data/astro_nearest_peak_to_gene_TSS.csv")

ns5.ids <- ns5[,"EnsemblID"]

astro.ids <- astro[,"EnsemblID"]

shared <- merge(ns5,astro, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_NS5","_Astro"))


##tidy up abit...

shared.tidy <- shared[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,34,35,36,37)]

write.csv(shared.tidy, file = "results/REST_peaks_shared_between_ns5_and_ns5dAstro.csv")

#ns5 unique peaks

ns5.unique.ids <- ns5.ids[which(!(ns5.ids %in% astro.ids))]

astro.unique.ids <- astro.ids[which(!(astro.ids %in% ns5.ids))]

#find data

ns5.unique <- ns5[which(ns5[,"EnsemblID"] %in% ns5.unique.ids),]

astro.unique <- astro[which(astro[,"EnsemblID"] %in% astro.unique.ids),]

ns5.unique.o <- ns5.unique[order(ns5.unique[,"neg10log10pVal"],decreasing = TRUE),]
astro.unique.o <- astro.unique[order(astro.unique[,"neg10log10pVal"],decreasing = TRUE),]


write.csv(ns5.unique.o, file = "results/REST_peaks_unique_to_ns5.csv")
write.csv(astro.unique.o, file = "results/REST_peaks_unique_to_astro.csv")


################link back to expression changes when KO REST with domneg - will need to treat with caution tho as Kee Yew originally filtered and background corrected this file - and its been annotated with REMOAT which is a wee bit dirty...

limma <- read.csv(file="data/limma_rd.csv")

#massively tidy this up...keep chr positions, FC, Pval,genbank transcripts (pull out name), and ensembl id

keep.cols <- qw(Chr, Start, End, Strand, logFC, P.Value, adj.P.Val,Proportion_UCSC_transcripts,Proportion_Ensembl_transcripts)

limma.tidy <- limma[,keep.cols]

#remove bits in brackets

ucsc <- as.character(limma.tidy[,"Proportion_UCSC_transcripts"])

#try  gsub again...

ucsc.brackets <- substr(gsub('.*\\(','',ucsc),0,nchar(gsub('.*\\(','',ucsc))-1)

ensembl <- as.character(limma.tidy[,"Proportion_Ensembl_transcripts"])

ensembl.brackets <- substr(gsub('.*\\(','',ensembl),0,nchar(gsub('.*\\(','',ensembl))-1)

new.limma <- cbind(limma.tidy[,c(1,2,3,4,5,6,7)], ucsc.brackets, ensembl.brackets)

#new column names

colnames(new.limma) <- qw(Chr, Start, End, Strand, logFC, P.Value, adj.P.Val, Symbol, EnsemblID)


##remove duplicates - doing everything on symbol IDs as quite a few genes in the limma annotation do not have ensembl IDs....
#actually doing it on EnsemblID as it gives more matches between the files...

new.limma.o <- new.limma[order(abs(new.limma[,"logFC"]),decreasing = TRUE),]
new.limma.od <- new.limma.o[!duplicated(new.limma.o[,"EnsemblID"]),]

#remove gene lines that don't have a symbol annotation (\\w regexp takes any alphanumeric character)

limma.res <- new.limma.od[grep("\\w",new.limma.od[,"EnsemblID"]),]

#leaves us with 10875 genes

##merge with astrocyte peak data

merged <- merge(limma.res, astro, by.x = "EnsemblID", by.y = "EnsemblID", all.x = TRUE, suffixes=c("_limma","_IP"))

###try to mark genes that are unique with a yes/no answer - astro.unique.ids

unique.astro.peaks <- ifelse(merged[,"EnsemblID"] %in% astro.unique.ids, "Yes", "No")

peakexp <- cbind(merged,unique.astro.peaks)

#tidy up abit

keep.cols <- c(1,9,2,3,4,5,6,7,8,11,15,16,17,18,19,23,24,25,26,28,29)

peakexp.results <- peakexp[,keep.cols]

#order and save
peakexp.results.o <- peakexp.results[order(peakexp.results[,"logFC"],decreasing = TRUE),]
write.csv(peakexp.results.o, file = "results/DNREST_in_NS5dAstros_with_REST_peaks.csv")

##In total, studying 10387 genes

# of which 743 map to a REST peak (728 are unique to Astrocytes)

# do something with gene expression

#cut off below adj.P.Val <= 0.05

results <- peakexp.results.o[which(peakexp.results.o[,"adj.P.Val"] <= 0.05),]

#cut off FC above +/- 1

results.change <- results[which(abs(results[,"logFC"]) >= 1),]

write.csv(results.change,file = "results/DNREST_in_NS5dAstros_with_REST_peaks_expression_changes.csv")

###genes that change that have REST bound:
# [1] Acsl6         Dcn           Lrp11         Hcn2          Adora2b
# [6] 1700020C11Rik Armc2         Tle6          Pib5pa        Tph2
#[11] 6330407J23Rik Lgr5          Agxt2l2       Irf1          Phyhipl
#[16] mKIAA1981     Sunc1         Lyrm7         Pcbp3         Ddt
#[21] Ggtla1        Dppa1         Txnrd1        Pemt          Aifm2
#[26] Ndg2          Fyn           Igfbp3        BC005764

#upregulated

up.res <- results[which(results[,"logFC"]>=1),]

write.csv(up.res, file = "results/DNREST_in_NS5dAstros_with_REST_peaks_upregulated.csv")

###genes that are upregulated that have REST bound:

# [1] Acsl6         Dcn           Lrp11         Hcn2          Adora2b
# [6] 1700020C11Rik Armc2         Tle6          Pib5pa        Tph2
#[11] 6330407J23Rik Lgr5          Agxt2l2       Irf1          Phyhipl
#[16] mKIAA1981     Sunc1         Lyrm7         Pcbp3         Ddt
#[21] Ggtla1

#downregulated

down.res <- results[which(results[,"logFC"]<=-1),]

write.csv(down.res, file = "results/DNREST_in_NS5dAstros_with_REST_peaks_downregulated.csv")

###genes that are downregulated that have REST bound:

#[1] Dppa1    Txnrd1   Pemt     Aifm2    Ndg2     Fyn      Igfbp3   BC005764

###look for enrichment in genes that change that have REST bound...

upregulated <- length(up.res[,"EnsemblID"])
downregulated <- length(down.res[,"EnsemblID"])

uppeaks <- length(which(!is.na(up.res[,"Peak"])))
downpeaks <- length(which(!is.na(down.res[,"Peak"])))

allgenes <- length(merged[,"EnsemblID"])
allpeaks <- length(which(!is.na(merged[,"Peak"])))

matrix <- matrix(c(uppeaks,upregulated,allpeaks,allgenes),nr = 2)

fisher.test(matrix)

upregulated genes are not sig enriched for rest peaks (pval 0.5455)
neither are the down regulated (pval 1)






























