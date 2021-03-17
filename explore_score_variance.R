#run_locally <- T
#source("recombine.R")
library(Hmisc)
source("helpers.R")
library(gridExtra)
library(glue)

CROHNS <- "crohns"

SCHIZ <- "schiz"
SCHIZ_PREVALENCE <- 0.01

MakePanelFigure <- function(crohns_prevalence = 0.005) {
  l1 <- variance_plot(CROHNS, crohns_prevalence, alpha = 0.04)
  l2 <- variance_plot(SCHIZ, SCHIZ_PREVALENCE, alpha = 0.04)
  pdf(glue("Results/variance_plot_crohns_prevalence{crohns_prevalence}.pdf"), height = 6, width = 6)
  grid.arrange(l2$p1 + ggtitle("A (Schizophrenia)"), l2$p2 + ggtitle("B (Schizophrenia)") + theme(legend.position = "none"), 
               l1$p1 + ggtitle("C (Crohn's  Disease)"), l1$p2 + ggtitle("D (Crohn's Disease)") + theme(legend.position = "bottom"), 
               nrow = 2)
  dev.off()
}

variance_plot <- function(disease, prevalence, alpha) {
  l <- get_data(disease = disease, prevalence = prevalence, num_couples = 5000)
  children <- l$sim_scores
  sampled <- children %>%
    group_by(p1, p2) %>%
    summarise(var = var(score), 
              mean = mean(score), 
              average_parental_score = first((p1_score + p2_score) / 2), 
              num_parents_cases = first((p1_pheno == 1) + (p2_pheno == 1)), 
              .groups = "drop")
  
  ee <- ecdf(sampled$average_parental_score)
  sampled$q <- ee(sampled$average_parental_score)
  
  
  parents <- l$p_info %>%
    group_by(pheno) %>%
    mutate(weight = if_else(pheno == 0, (1 - prevalence) / n(), prevalence / n()))
  wtd_var <- wtd.var(parents$score, weights = parents$weight, normwt = T)
  print(disease)
  cat("Variance parents (weighted by phenotype):", wtd_var, "\n")
  #cat("Unweighted variance parents:", var(parents$score), "\n")
  cat("Variance children (sampled based on parental phenotype):", var(children$score), "\n")
  cat("Mean within-family variance (after sampling based on parental phenotype):", mean(sampled$var), "\n")
  
  
  p1 <- ggplot(sampled, aes(x = q, y = var)) + 
    geom_smooth() +#method = 'gam',
    #formula = y ~ splines::ns(x, 5)) + #ylim(0, 5e-09) + 
    xlab("Mid-parental PRS (quantile)") + 
    ylab("Variance of embryo PRSs") + theme_bw(13) + 
    geom_point(alpha = alpha)+ 
    geom_hline(yintercept = wtd_var, alpha = 0.25) + 
    coord_cartesian(xlim = c(0, 1), ylim=c(0, wtd_var * 1.5), expand = F) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", 1))
    
  n <- 20
  sims <- rchisq(600000, df = n - 1) / (n-1) * (wtd_var / 2)
  
  for_plot <- bind_rows(sampled %>% transmute(type = "observed", var), 
                  data.frame(var = sims) %>%
                    mutate(type = "theoretical"))
  p2 <- ggplot(for_plot, aes(x = var, col = type)) + 
    geom_density() + 
    xlab("Variance of embryo PRSs") +
    #geom_hline(yintercept = wtd_var) + 
    theme_bw(13) + 
    xlim(0, wtd_var * 2) + 
    ylab("Density") + 
    scale_color_discrete("")
  
  return(list(p1 = p1, p2 = p2))
  
}

# for (id in sample_ids[1:500]) {
#   
#   parent_vcf <- GetParentVCF(vcf = vcf, sample_id = id)
#   
#   score <- mean(0.5 * (1 - as.numeric(merged$h1)) * merged$BETA + 0.5 * (1 - as.numeric(merged$h2)) * merged$BETA)
#   
#   merged[, `:=`(het = !h1 == h2)]
#   merged[, `:=`(var = if_else(het, 0.25 * 0.25 * merged$BETA ^ 2, 0))]
#   total_var <- sum(merged$var) / nrow(merged)^2
#   vars <- c(vars, total_var)
# }


get_var_contributions <- function() {
  if (F) {
    vcf_o <- fread(VCF_FILE, skip = "#CHROM")
    vcf <- copy(vcf_o)
    setnames(vcf, '#CHROM', "CHROM")
  }
  gwas <- fread(GWAS_FILE)[, .(ID, new_BETA, REF, ALT)]
  setnames(gwas, "new_BETA", "BETA")
  freqs <- fread(FREQ_FILE)
  
  joined <- merge(gwas, freqs, by = "ID")
  ggplot(joined, aes(x = 1 - ALT_FREQS, y = BETA)) + 
    geom_point() +
    xlab("Frequency effect allele") + 
    theme_bw()
  ggsave("beta_versus_freq.pdf")
  
  ggplot(joined, aes(x = BETA)) +
    geom_histogram(bins=1000) + 
    theme_bw() + 
    xlab("Beta")
  ggsave("beta_histogram.pdf")
  
  joined$abs_beta <- abs(joined$BETA)
  joined$beta_neg <- joined$BETA < 0
  joined$effect_allele_freq <- 1 - joined$ALT_FREQS
  
  sample_ids <- colnames(vcf)[10:ncol(vcf)]
  dir <- "~/Documents/Docs/embryos/"
  scores_file <- glue("{dir}/LIJMC_scores.sscore")
  
  scores <- fread(scores_file)
  setnames(scores, '#IID', 'ID')
  other_dat <- GetAges()[, .(ID, plink_pheno)]
  sel_cols <- vcf[, c("ID", sample_ids), with = F]
  melted <- melt(sel_cols, id.vars = "ID", measure.vars = sample_ids)
  sep <- melted[value %in% c("0|1", "1|0"), ]
  sep[, c("h1", "h2") := tstrsplit(value, "|", fixed = T)]
  merged <- merge(sep, gwas, by = "ID")
  summ <- merged[, .(n = .N, var_contrib = sum(BETA^2)), by = variable]
  setnames(summ, 'variable', "ID")
  var_contrib0 <- merge(summ, scores, by = "ID")
  var_contrib <- merge(var_contrib0, other_dat, by = "ID")
  return(var_contrib)
  ggplot(mmm, aes(x = SCORE1_AVG, y=var_contrib)) + 
    geom_smooth() + 
    facet_wrap(~plink_pheno) + 
    xlab("PRS") + 
    geom_point() +
    ylab("PRS variance contribution")
  ggsave("variance_contribution_by_prs.pdf")
  ggplot(mmm, aes(x = SCORE1_AVG, y=n)) + 
    geom_smooth() + 
    facet_wrap(~plink_pheno) + 
    xlab("PRS") + 
    geom_point() +
    ylab("Num het SNPs")
  ggsave("num_het_snps_by_prs.pdf")
  
  var_contrib <- map(sample_ids[1:1000], function(s){
    vcf <- GetParentVCF(vcf = vcf, sample_id = s)
    merged <- merge(vcf, gwas, by = "ID")
    hets <- merged[!h1==h2, ]
    var_contrib <- sum(hets$BETA^2)
    data.frame(`#IID` = s, var_contrib = var_contrib, num_hets = nrow(hets))
  })
  browser()
  
  vars <- c()
  map2 <- fread(glue("{MAP_DIR}/sexavg_chr2.txt"))
  id <- 'GCOBBN4043'
  parent_vcf <- GetParentVCF(vcf = vcf, sample_id = id)
  
  grid <- expand.grid(factor = seq(0.5, 5, by = 0.5), 
                      num_kid = 1:10000)
  
  parent_chr <- parent_vcf[CHROM == 2, ]
  setorder(parent_chr, POS)
  merged <- merge(parent_chr, gwas, by = 'ID')
  setorder(merged, POS)
  browser()
  res <- pmap_dfr(grid, function(factor, num_kid){
    map <- map2 %>%
      mutate(cM = cM * factor)
    which_haplotype <- GetHaplotypeAssignment(map, parent_chr$POS)
    passed <- as.numeric(GetHaplotype(parent_chr, which_haplotype))
    score <- mean(passed * merged$BETA)
    data.frame(factor = factor, num_kid = num_kid, score = score)
  })
  summ <- group_by(res, factor) %>% summarise(s = var(score))
  p <- ggplot(summ, aes(x=factor, y=s)) + 
    geom_line() + 
    xlab("Factor multiplying genetic distances") + 
    ylab("Variance of gamete PRS")
  browser()
}
