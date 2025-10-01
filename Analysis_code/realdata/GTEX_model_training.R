library(SMiXcan)

# SMiXcan ------------------------------------------------------------------


setwd("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data")
# load GTEx covariate data
cov1=data.frame(fread("phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"))
cov0=cov1[,c("SUBJID", "AGE")]
cov2=fread("Breast_Mammary_Tissue.v8.covariates.txt")
cov3=t(cov2[,-1])
colnames(cov3)=data.frame(cov2)[,1]
cov3=data.frame(colnames(cov2)[-1], cov3);colnames(cov3)[1]="SUBJID"
cov3=cov3[,c(1:21,67:69)] # top 15 PEER factors
cov=data.frame(merge(cov0, cov3, by="SUBJID"))
cov <- cov%>%
  filter(SUBJID %in% gtex_white) %>%
  filter(sex==2) %>%
  select(!sex)

#> dim(cov)
#[1] 125  24


# load genotype
geno1=fread("shapeit_data_for_predictdb_variants-r2") 
# 178698    847
geno=geno1[,c(1:9,match(cov[,1], colnames(geno1))), with=F]
#> dim(geno)
#[1] 178698    134
#
# Load model: elastic net model output
filename <- "en_Breast_Mammary_Tissue.db"
library(DBI)
sqlite.driver <- dbDriver("SQLite")
ElasticNet <- dbConnect(sqlite.driver, dbname = filename)
dbListTables(ElasticNet)
ENextra <- dbReadTable(ElasticNet,"extra")
ENweights <- dbReadTable(ElasticNet,"weights")

# overlapping genes
genID=genID1=intersect(ENextra$gene, breast$gene_id)
length(genID)
# 6461 correct


#run (S)-MiXcan GTEx training - fix seed.----------------

result <- vector("list", 6461)
set.seed(13313)
t1=Sys.time()


# results container
res_weights_all <- vector("list", length = 6461)

# main loop

for (j in 1:6461){
  print(j)
  # indicator
  yName=genID[j]
  xName=ENweights[which(ENweights$gene==yName), "varID"]
  xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
  nName=cov$"SUBJID" # women
  n=length(nName)
  # TODO: clean below ( Extract the 125 European  GTEx sample for gene j
  yData=t(breast[which(breast$gene_id==yName), match(nName, colnames(breast)), with=F])
  xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with=F])
  zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
  piData=pis[match(nName, rownames(pis)),2]
  xCellscore.epiData=xCellscore.epi[match(nName, xCellscore.epi[,1]), 2]
  class(xData)<-"numeric"
  cp.idx=complete.cases(xData) & complete.cases(yData)
  
  px=ncol(xData)
  if (px>1) {
    xvar0=which(apply(xData[cp.idx,], 2, function(f) mean(f)>0.05))
    x.complete=xData[cp.idx,xvar0]
    #xz.test=xz.test1[,c(xvar0, (px+1):ncol(xz.complete1))]
  }
  if (px==1) {
    xvar0=1*(mean(xData[cp.idx,])>0.05)
    x.complete=matrix(xData[,xvar0])
  }
  if (ncol(x.complete) == 0 ||is.null(nrow(x.complete))) {next}
  
  px2=ncol(x.complete)
  z.complete=zData[cp.idx,]
  xz.complete=as.matrix(cbind(x.complete, z.complete))
  y.complete=yData[cp.idx]
  pi.complete=piData[cp.idx]
  length(y.complete)

  pz=ncol(z.complete)
  # 10 fold cross validation
  set.seed(1334 + j*149053)
  foldid= sample(1:10, length(y.complete), replace=T)
  
  #   # MiXcan method
  # to fix training_prediction_model in S-MiXcan
  ft.sym=SMiXcan::MiXcan(y=y.complete, x=x.complete, cov=z.complete, pi=pi.complete, xNameMatrix=xName.all[xvar0,], foldid=foldid)
 
    # combine weights for both cells on the same SNP rows
    w1 <- ft.sym$beta.SNP.cell1 %>% rename(weight_cell_1 = weight)
    w2 <- ft.sym$beta.SNP.cell2 %>% rename(weight_cell_2 = weight)
    
    w <- inner_join(
      w1, w2,
      by = c("gene","rsid","varID","ref_allele","eff_allele")
    ) %>%
      mutate(
        gene_id   = gene_id_stripped,
        gene_name = gene_name,
        type      = ft.sym$type
      ) %>%
      select(
        gene_id, gene_name, gene, rsid, varID, ref_allele, eff_allele,
        weight_cell_1, weight_cell_2, type
      )
    
    if (nrow(w)) res_weights_all[[j]] <- w
    if (j %% 200 == 0) cat("Processed", j, "genes\n")
  }
  
  #  bind & save after loop
  weights_final <- bind_rows(res_weights_all)
  
  
  #Save result
  write_csv(weights_final, "weights_miXcan_full_2025.csv")

  print(dim(weights_final))
  head(weights_final, 10)
  
  