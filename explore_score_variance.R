run_locally <- T
source("recombine.R")
source("helpers.R")
source("get_ages.R")

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


compare_variance <- function(remove_dups = F) {
  l <- get_data_helper()
  children <- l$sim_scores
  var_children <- children %>%
    group_by(p1, p2) %>%
    summarise(var = var(score), 
              mean = mean(score))
  p_info <- l$p_info
  setnames(p_info, c("#IID", "pheno"))
  parents <- l$parents_scores
  parents <- merge(parents, p_info, by = "#IID")
  
  parent_scores0 <- unique(children[, .(p1, p2)])
  parents1 <- parents %>%
    transmute(p1 = `#IID`, p1_score = SCORE1_AVG, p1_pheno = pheno)
  parents2 <- parents %>%
    transmute(p2 = `#IID`, p2_score = SCORE1_AVG, p2_pheno = pheno)
  parent_scores <- reduce(list(parent_scores0, parents1, parents2), left_join) %>%
    mutate(average_score = (p1_score + p2_score) / 2) 
  bins_data <- left_join(var_children, parent_scores, by = c("p1", "p2")) %>%
    mutate(num_parents_schiz = (p1_pheno == 2) + (p2_pheno == 2), 
           parent_score_diff = abs(p1_score - p2_score)) %>%
    ungroup()
  # variance depends on difference between parent scores as well as average score
  # difference depends on average score
  mod <- lm(var ~ poly(parent_score_diff, 3), bins_data)
  bd <- bins_data %>%
    ungroup() %>%
    mutate(resid = residuals(mod))
  ggplot(bd, aes(x = average_score, y = resid)) + geom_smooth()
  
  if (remove_dups) {
    bins_data <- group_by(bins_data, p1) %>%
      slice(1) %>%
      group_by(p2) %>%
      slice(1)
  }
  
  bins_data <- bins_data %>%
    group_by(num_parents_schiz) %>%
    mutate(max = case_when(num_parents_schiz == 0 ~ if_else(remove_dups, 962, 10000),   #962
                           num_parents_schiz == 1 ~ if_else(remove_dups, 200, 19), #200, # 19,      
                           num_parents_schiz == 2 ~ 1)) #0))      # 1
  
  sampled <- bins_data %>%
    group_by(num_parents_schiz) %>%
    mutate(i = 1:n()) %>%
    filter(i <= max)
  ee <- ecdf(sampled$average_score)
  sampled$q <- ee(sampled$average_score)
  p <- ggplot(sampled, aes(x = q, y = var)) + 
    geom_smooth() + #ylim(0, 5e-09) + 
    xlab("Mid-parental PRS (quantile)") + 
    ylab("Variance of embryos PRSs") + theme_bw() + geom_point(alpha = 0.1)#+ 
  #facet_wrap(~num_parents_schiz)
  print(p)
  ggsave("variance_by_midparental_prs_quantiles.pdf")
  browser()
  
  parents <- parents %>%
    group_by(pheno) %>%
    mutate(weight = if_else(pheno == 1, 0.99 / n(), 0.01 / n()))
  wtd_var <- wtd.var(parents$SCORE1_AVG, weights = parents$weight, normwt = T)
  cat("Weighted variance parents:", wtd_var, "\n")
  cat("Unweighted variance parents:", var(parents$SCORE1_AVG), "\n")
  cat("Mean variance children conditional on parents:", mean(var_children$var), "\n")
  cat("Variance parents by phenotype\n")
  parents %>%
    ungroup() %>%
    mutate(pheno = if_else(pheno==2, "Schiz.", "Not schiz.")) %>%
    group_by(pheno) %>%
    summarise(v = var(SCORE1_AVG))
  
  var <- wtd_var
  n <- 20
  sims <- rchisq(600000, df = n - 1) / (n-1) * (wtd_var / 2)
  var_children$v2 <-  var_children$var
  aa <- bind_rows(var_children %>% transmute(type = "observed", var = v2), 
                  data.frame(var = sims) %>%
                    mutate(type = "theoretical"))
  ggplot(aa, aes(x = var, col = type)) + 
    geom_density() + 
    xlab("Sampling variance") +
    geom_vline(xintercept = wtd_var / 2)
  
  browser()
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

