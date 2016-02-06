## Copyright 2016 Arjun Arkal Rao
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
args <- commandArgs(trailingOnly = TRUE)
library(tools)

print_usage = function(){
    print(paste("Usage: Rscript mhci_rankboost.R <MHC_pred_file> <rsem_isoforms_file> ",
              "<transgened_peptide_fasta> [<V,W,X,Y,Z>]", sep="", collapse=""))
    print(paste("Where V,W,X,Y,Z = boost values for num_pepts_all, num_pepts_high, num_mhc, tpm ",
                "respectively. The values should be comma separated without any spaces.",
                sep="", collapse=""))
    quit("no", 1) 
}

if (length(args) == 3){
        ratios=c(1, 1, 1, 1)
} else if (length(args) == 4){
    ratios <- as.numeric(strsplit(args[4], ",")[[1]])
    if (length(ratios) != 4 ) {
        print("ERROR: Need to submit 4 values for ratio") 
        print_usage()
    }
} else {
    print_usage()
}

filename <- args[1] # input merged mhc predicitons file
rsem_fn <- args[2] # input rsem call containing file for the sample
fa_fn <- args[3] # input peptide file for the prediction

max_boost <- 55  # This allows peptide #2 to overthrow peptide #1
boost_npa <- max_boost * ( ratios[1]/sum(ratios) )
boost_nph <- max_boost * ( ratios[2]/sum(ratios) )
boost_nMHC <- max_boost * ( ratios[3]/sum(ratios) )
boost_TPM <- max_boost * ( ratios[4]/sum(ratios) ) 

# read in merged MHC file
x <- read.table(filename, colClasses=c("character", "character", "character", "NULL", "NULL",
                                       "numeric", "character", "character", "character"))
colnames(x) <- c("allele", "MHC_peptide", "input_peptide", "pc_binder", "ensembl_gene", "HUGO_gene",
                 "mutations")

#read in peptide file
fa_file<-read.table(fa_fn, header=F, colClasses=c("character"))
fa_file[, 1]<-gsub(">", "", fa_file[, 1])
# read in the isoform level expression file.  Then split it by gene name
rsem <- read.table(rsem_fn, header=T, row.names=1)

# Making a change here for easy printing later
x <- x[, c("allele", "input_peptide", "MHC_peptide", "pc_binder", "mutations", "ensembl_gene",
           "HUGO_gene")]

#split by input IAR
all_data <- split(x, x$input_peptide)
#get min PC for each IAR
stats <- as.data.frame(sapply(all_data, function(x){min(x$pc_binder)}))
colnames(stats)[1] <- "min_pc"
#number of peptides within this IAR
stats$"num_pept" <- sapply(all_data, function(x){length(x[, 1])})
#number of MHCs binding to peptides with this core
stats$"num_MHC" <- sapply(all_data, function(x){length(split(x, x$allele))})

stats$gene <- sapply(all_data, function(x){x[1,"HUGO_gene"]})
stats$ensgene <- sapply(all_data, function(x){x[1,"ensembl_gene"]})

stats$mutation <- sapply(all_data, function(x){
    z <- strsplit(as.character(x[1, "mutations"]), ",")[[1]]
    z <- strsplit(z, '_')[[1]][-1]
    mutations <- gsub('[1234567890]+', '>', z)
    positions <- as.numeric(gsub('[A-Z]+', '', z))
    outstring <- mutations[1]
    if (length(mutations) > 1){
        for(i in c(2:length(mutations))){
            outstring <- paste(outstring, "+", positions[i]-positions[i-1]-1,
                               "X+", mutations[i], sep='', collapse='')
        }}
    outstring
})

stats$peptide <- sapply(rownames(stats), function(x){
    fa_file[grep(paste('^', x, '$', sep='', collapse=''), fa_file[,1])+1, 1]})

# for x, find the values from val_col (5 for TPM, 6 for FPKM) for the expression
get_rsem_val <- function(stats, rsem_col){
    sapply(rownames(stats), function(x){
        z <- strsplit(as.character(all_data[x][[1]][1,"mutations"]), ",")[[1]]
        z <- gsub("_.*", "", z)
        sum(as.numeric(rsem[z,rsem_col]))
    })}

stats$TPM <- get_rsem_val(stats, 5)

stats<-stats[stats$TPM!=0&!is.na(stats$TPM),]
if ( length(stats[,1]) ==0 ){
  print("No mutations found in expressed genes")
  quit("no", 1)
}
medTPM <- median(rsem[rsem[, 5]!=0,5])
stats <- stats[stats$TPM > (0.1 * medTPM), ] # this can change
if ( length(stats[,1]) ==0 ){
  print("No mutations found in expressed genes")
  quit("no", 1)
}

stats<-stats[order(stats$min_pc, - stats$num_pept/stats$num_MHC, -stats$TPM),]

stats$old_rank <- c(1: length(stats[, 1]))
stats$mod_rank <- c(1: length(stats[, 1]))

med_MHC<-median(stats$num_MHC)
MHC_boost<- max(stats$num_MHC) - c(max(stats$num_MHC):median(stats$num_MHC))[-1]
print("stats successfully created")

if (length(stats[,1]) > 10) {
  for (i in rownames(stats)[-1]) {
    neoepitope <- rownames(stats[i,])
    # boost by num_pept_high
    n <- sum(all_data[neoepitope][[1]][,"pc_binder"] <= (stats[i, 'min_pc'] + 0.3)) - 1 
    boost <- boost_nph * ( (n >= 1) * 0.4 + 
                            (n >= 2) * 0.3 + 
                            (n >= 3) * 0.2 + 
                            (n >= 4) * 0.1 )
    # boost by num_pept_all
    n <- length(all_data[neoepitope][[1]][,"pc_binder"])-1 
    boost <- boost_npa * ( (n >= 10) * 0.4 +  
                            (n >= 15) * 0.3 +  
                            (n >= 20) * 0.2 +  
                            (n >= 30) * 0.1 ) 
    # boost by nMHC
    boost_level <- stats[i, "num_MHC"] - med_MHC
    boost <- boost + boost_nMHC * 
             ifelse( boost_level > 0, sum(MHC_boost[c(1:boost_level)]) / sum(MHC_boost), 0 ) 
    # boost by TPM
    boost <- boost + boost_TPM * ( (stats[i, "TPM"] > medTPM) * 0.35 + 
                                   (stats[i, "TPM"] > (2 * medTPM)) * 0.25 + 
                                   (stats[i, "TPM"]> (5 * medTPM)) * 0.15 +
                                   (stats[i, "TPM"]> (10 * medTPM)) * 0.25 )

    oldrank <- stats[i, 'mod_rank']
    newrank <- round((100 - boost) * stats[i, 'mod_rank'] / 100)
    if(oldrank == newrank){
      next
    }
    stats <- rbind(
        stats[stats$mod_rank<newrank, ],
        stats[stats$mod_rank==oldrank, ],
        stats[stats$mod_rank >= newrank & stats$mod_rank < oldrank, ],
        stats[stats$mod_rank>oldrank, ]
        )
    stats$mod_rank <- c(1: length(stats[, 1])) 
  }
  stats<-stats[order(stats$mod_rank),]
  stats<-stats[stats$TPM>quantile(stats$TPM,0.25),]
  stats$mod_rank<-c(1:length(stats$TPM))
  print("loop completed")
} else{
  print("Too few calls to run the loop. Leaving in original sorted mode.")
}

stats$"binding_MHCs" <- sapply(rownames(stats), function(x){
  paste(unique(unlist(sapply(all_data[x], "[[",1))),sep=",", collapse=",")
})

file_prefix <- gsub(paste(".", file_ext(filename), sep="", collapse=""), "", filename)
write.table(paste("gene", "mutation", "containing_peptide", "corrected_rank", "original_rank",
                  "num_binding_MHC", "filter_passing_peptides", "min_percentile_binder", "TPM",
                  "binding_MHCS" , sep="\t", collapse="\t"),
            file=paste(file_prefix, "_concise_results.tsv", sep="", collapse=""),
            col.names=F, row.names=F, quote=F)

write.table(file=paste(file_prefix, "_concise_results.tsv", sep="", collapse=""),
            stats[,c("gene", "mutation", "peptide", "mod_rank", "old_rank", "num_MHC", "num_pept",
                     "min_pc", "TPM", "binding_MHCs")], col.names=F, row.names=F,
            append=T, quote=F, sep="\t")

for (neoepitope in rownames(stats))
{
  write.table(file=paste(file_prefix,"_detailed_results.tsv", sep="", collapse=""),
              x=paste("#", stats[neoepitope, "gene"], "\t", stats[neoepitope,"mutation"], "\t",
                      stats[neoepitope, "peptide"], sep="", collapse=""), col.names=F, row.names=F,
              quote=F, append=T)
  columns <- colnames(all_data[neoepitope][[1]])
  columns <- columns[!columns%in%c("HUGO_gene", "ensembl_gene")]
  write.table(paste(columns, sep="\t", collapse="\t"), row.names=F, col.names=F, quote=F, append=T,
              sep="\t", file=paste(file_prefix,"_detailed_results.tsv", sep="", collapse=""))
  write.table(all_data[neoepitope][[1]][order(all_data[neoepitope][[1]][,"pc_binder"]),columns],
              row.names=F, col.names=F, quote=F, append=T, sep="\t",
              file=paste(file_prefix, "_detailed_results.tsv", sep="", collapse=""))
  write.table(file=paste(file_prefix, "_detailed_results.tsv", sep="", collapse=""), x="\n",
              col.names=F, row.names=F, quote=F, append=T)
}

