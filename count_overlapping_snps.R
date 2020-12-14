library(data.table)
library(glue)
library(purrr)

daner <- fread("gunzip -c /vol/sci/bio/data/shai.carmi/db2175/embryo_selection/daner_PGC_SCZ_w3_90_0518d__lencz.gz")
dir <- "/vol/sci/bio/data/shai.carmi/db2175/embryo_selection/LIJMC/"

res <- map_dfr(1:22, 
  function(x){
		daner_chr <- daner[CHR == x, ]
		num_daner_snps <- nrow(daner_chr)
		map_chr <- fread(glue("{dir}LIJMC37_.phased.cm.{x}.map"))
		num_lijmc_snps <- nrow(map_chr)
		daner_map_chr <- daner_chr[SNP %in% map_chr$V2, ]
		num_common_snps <- nrow(daner_map_chr)
		daner_map_chr_0.05 <- daner_map_chr[P <= 0.05, ]
		num_snps_lt_0.05 <- nrow(daner_map_chr_0.05)
		clumps <- fread(glue("{dir}/clumped/{x}.clumped"))
		num_clumps <- nrow(clumps)
		data.frame(chr=x, num_daner_snps = num_daner_snps, 
		           num_lijmc_snps = num_lijmc_snps, 
							 num_common_snps = num_common_snps, 
							 num_snps_lt_0.05 = num_snps_lt_0.05, 
							 num_clumps = num_clumps)
	}
	)
fwrite(res, "clump_qc.csv")
