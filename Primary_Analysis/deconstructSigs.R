library(deconstructSigs)
library(RColorBrewer)

par(mar=c(1,1,1,1))
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
items_to_plot <- args[2]


workDir <- file.path(path,"SigMutations")
#workDir <- "/Users/luol2/Dropbox/MSK/Projects/Ben_Project/WES_mouse_Project_10212_F/SigMutations"
dir.create(workDir, showWarnings = FALSE)

#items <- read_tsv(file = file.path(path, "items_to_plots_MutSig.txt"))
items <- read.table(file = file.path(path, items_to_plot), sep = "\t", header = T, quote = "")

tumors <- apply(items[,1:2], 1, paste, collapse = ".")


for (i in 1:length(tumors)) {
    t <- tumors[i]
    vcf <- file.path(path, items$files[i])
    #vcf <- list.files(path = paste0(workDir,'/',t), pattern = "*.vcf", full.names =  TRUE)
    sigs.input <- vcf.to.sigs.input(vcf = vcf)
    wchSig <- whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.exome.cosmic.v3.may2019, contexts.needed = TRUE, tri.counts.method = 'exome')
    #pdf(file = file.path(workDir, paste0(t,'.pdf')),  width = 10, height = 10)
    png(file = file.path(workDir, paste0(t,'.png')), height = 800, width = 800)
    #chart <- plotSignatures(wchSig, sub = t)
    chart <- plotSignatures(wchSig)
    dev.off()

    n <- length(which(wchSig$weights>0))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    #pdf(file = file.path(workDir, paste0(t,'.piechart.pdf')), width = 10, height = 10)
    png(file = file.path(workDir, paste0(t,'.piechart.png')), height = 800, width = 800)
    #pie <- makePie(wchSig, add.color = sample(col_vector,n), sub = t)
    pie <- makePie(wchSig, add.color = sample(col_vector,n))
    dev.off()

}
