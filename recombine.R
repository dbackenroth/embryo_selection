library(data.table)
library(glue)
library(dplyr)
library(purrr)
library(tidyr)

options(stringsAsFactors = F)

MAP_DIR <- "/Users/dbackenr/OneDrive - JNJ/huji/haploseek_v2/GeneticMap/Refined_genetic_map_b37/"
MAP_DIR <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/Maps/Refined_genetic_map_b37/"
OUT_DIR <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/Peds/"
COUPLES_FILE <- "//cs/icore/db2175/embryo_selection/selected_couples.csv"
VCF_FILE <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC_score_snps.recode.vcf"
CALCULATE_SCORES_SCRIPT <- "./calculate_scores_children.sh"

#couples <- data.frame(p1 = c("SZPABR0002", "SZPABR0003"), 
#                      p2 = c("SZPABR0003", "SZPABR0005"))

GetHaplotype <- function(chr_data, which_haplotype) {
  flip <- sample(c(T, F), 1)
  haplotypes <- as.matrix(chr_data[, .(h1, h2)])
  if (flip) {
    which_haplotype <- 3 - which_haplotype
  }
  indices <- cbind(1:nrow(haplotypes), 
                   which_haplotype)
  passed_haplotype <- haplotypes[indices]
}

GetParentVCF <- function(vcf, sample_id) {
  vcf_sel <- copy(vcf)
  setnames(vcf_sel, sample_id, "h")
  vcf_sel <- vcf_sel[, .(CHROM, POS, ID, h)]
  vcf_sel[, c("h1", "h2") := tstrsplit(h, "|", fixed = T)]
  return(vcf_sel)
}

GetHaplotypeAssignment <- function(map, pos) {
  interp_cm <- approx(x = map$pos, y = map$cM, xout = pos)$y
  if (is.na(interp_cm[1])) {
    first_non_na <- min(which(!is.na(interp_cm)))
    interp_cm[1:(first_non_na - 1)] <- min(interp_cm, na.rm = T)
  }
  if (is.na(interp_cm[length(interp_cm)])) {
    interp_cm[is.na(interp_cm)] <- max(interp_cm, na.rm = T)
  }
  
  min_cM <- min(interp_cm)
  max_cM <- max(interp_cm)
  poisson_mean <- (max_cM - min_cM) / 100
  n_recombs <- rpois(1, poisson_mean)
  recombs_cM <- sort(runif(n_recombs, min = min_cM, max = max_cM), 
                     decreasing = T)
  haplotype <- rep(1, length = length(interp_cm))
  which_haplotypes <- 1:(n_recombs) %% 2 + 1
  if (n_recombs > 0) {
    for (i in 1:length(which_haplotypes)) {
      haplotype[interp_cm <= recombs_cM[i]] <- which_haplotypes[i]
    }
  }
  r <- rle(haplotype)
  # if two recombination events are very close to each other then 
  # this won't be true
  #if (!length(r$lengths) == n_recombs + 1){
  #  browser()
  #}
  return(haplotype)
}

GenerateChildren <- function(grid, vcf) {
  
  children <- pmap_dfr(grid, function(chrom, child, p1, p2, vcf) {
    map <- fread(glue("{MAP_DIR}/sexavg_chr{chrom}.txt"))

    vcf_chr <- vcf[CHROM == chrom, ]
    
    father_chr <- GetParentVCF(vcf_chr, p1)
    mother_chr <- GetParentVCF(vcf_chr, p2)
    
    which_father_haplotype <- GetHaplotypeAssignment(map, father_chr$POS)
    which_mother_haplotype <- GetHaplotypeAssignment(map, father_chr$POS)
    
    father_passed <- GetHaplotype(father_chr, which_father_haplotype)
    mother_passed <- GetHaplotype(mother_chr, which_mother_haplotype)
    
    data.frame(genotype = paste(father_passed, mother_passed, sep = "|"), 
               pos = father_chr$POS,
               chr = chrom, 
               ID = father_chr$ID,
               sample = paste(p1, p2, child, sep = "_"))
    
  }, vcf = vcf)
  #CHROM     POS       ID REF ALT QUAL FILTER INFO FORMAT   h h1 h2
  #1:    19 3579000 GA015527   A   G    .   PASS    .     GT 1|0  1  0
  #2:    19 3579000   rs1443   T   C    .   PASS    .     GT 1|0  1  0
  children_spread <- spread(children, sample, genotype)
  return(children_spread)
}

MakeChildrenVCF <- function(couples, n_children = 2, chrs = 1:22, 
                            out_file = "temp.vcf") {
  
  vcf_orig <- fread(VCF_FILE, skip = "#CHROM")
	#vcf_orig <- fread(cmd = "gunzip -c ~/Documents/Docs/embryos/A.vcf.gz", skip = "#CHROM")
  #header <- fread(cmd = "gunzip -c ~/Documents/Docs/embryos/A.txt.gz", nrows = 30)
  vcf <- copy(vcf_orig)
  setnames(vcf, "#CHROM", "CHROM")
  
  grid <- merge(data.frame(chrom  = chrs), 
                data.frame(child = 1:n_children)) %>%
    merge(couples)
  
  children <- GenerateChildren(grid = grid, vcf = vcf)
  
  vcf_orig_sel <- vcf_orig[, 1:9, with = F]
  children_dt <- as.data.table(children)
  setnames(children_dt, c("pos", "chr"), c("POS", "#CHROM"))
  both <- merge(vcf_orig_sel, children_dt, sort = F)
  writeLines("##fileformat=VCFv4.1", out_file)
  write.table(both, out_file, sep = "\t", append = T, row.names = F, quote = F)
}

DoRun <- function(groups, number, all = F) {
  set.seed(number)
	couples_grid <- read.csv(COUPLES_FILE, stringsAsFactors = F)
  if (!all) {
    couples_grid <- couples_grid %>%
      mutate(i = 1:n()) %>%
      filter(i %% groups == number - 1) %>%
      dplyr::select(-i)
  }
  grid <- couples_grid[, c("ID.1", "ID.2")]
	colnames(grid) <- c("p1", "p2")
	print(grid)
	out_vcf <- glue("{OUT_DIR}/{groups}_{number}.vcf")
	MakeChildrenVCF(grid, n_children = 20, chrs = 1:22, 
                  out_file = out_vcf)
	system(glue("{CALCULATE_SCORES_SCRIPT} {out_vcf}"))
}

# ##fileformat=VCFv4.1
