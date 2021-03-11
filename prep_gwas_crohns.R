library(data.table)
library(glue)
dir <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/crohns/"
daner_file <- glue("{dir}/EUR.CD.gwas_info03_filtered.assoc")
daner <- fread(daner_file)
vcf_file <- glue("gunzip -c {dir}/phased_no_ambiguous.vcf.gz | cut -f1,2,3,4,5")
vcf <- fread(cmd = vcf_file, skip = "#CHROM")
clumped <- fread(glue("{dir}/EUR.CD_no_ambiguous.clumped"))
daner_clumped <- daner[SNP %in% clumped$SNP, ]
vcf_clumped <- vcf[ID %in% clumped$SNP, ]
daner_clumped_sel <- daner_clumped[, .(CHR, SNP, BP, A1, A2, P, OR, SE)]
daner_clumped_sel[, `:=`(ID = SNP, BETA=log(OR))]
merged <- merge(daner_clumped_sel, vcf_clumped)
compl <- c(A = 'T', C = 'G', G = 'C', T = 'A')
merged[, `:=`(flip = A1 == ALT | A1 == compl[ALT]), ]
merged[, `:=`(new_BETA = ifelse(!flip, BETA, -BETA))]
new_daner <- merged[, .(ID, CHR, SNP, BP, REF, ALT, P, SE, new_BETA, flip)]
new_daner_file <- glue("{dir}/daner_flipped.tsv")
fwrite(new_daner, new_daner_file, sep = "\t") 
