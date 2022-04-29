#After LD calculation, PCA (Q) matrix, and relationship (K) matrix:

# get phenotype file and determine phenotype of interest
```{r, eval=FALSE}
pheno <- read.delim("http://ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", header=TRUE)
write.table(pheno[, -1:-2], "pheno.txt", 
            sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
```


---
  
  # GWAS using the gemma software package

```{bash, eval=FALSE}
#for other trait, change = 12(column for plant height in file), trait name in this pheno file
#Protein.content = column 36 in original file, 34 in pheno.txt
gemma -bfile binary_sativas413 -c pc3.txt -k output/binary_sativas413.cXX.txt -p pheno.txt -lmm 4 -n 34 -o Protein.content -miss 0.9 -r2 1 -hwe 0 -maf 0.05

cp output/Protein.content.assoc.txt cache
```

- `lmm`: specify frequentist analysis choice (default 1; valid value 1-4; 1: Wald test; 2:
                                                likelihood ratio test; 3: score test; 4: all 1-3.)
- `n`: specify phenotype column in the phenotype file (default 1); or to specify which
phenotypes are used in the mvLMM analysis
- `o`: specify output file prefix
- `miss`: specify missingness threshold (default 0.05)
- `r2`: specify r-squared threshold (default 0.9999)
- `hwe`: specify HWE test p value threshold (default 0; no test)
- `maf`: specify minor allele frequency threshold (default 0.01)

---
  
  # The Manhattan plot
  
  
  ```{r, eval=FALSE}
library(qqman)
library("data.table")
res <- fread("cache/Protein.content.assoc.txt")

pdf("graphs/Proteinmanhattan.pdf", width=10, height=10)
manhattan(x = res, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)
dev.off()
