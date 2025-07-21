## Analysis of Mintflow outputs for renal cell carcinoma dataset
## Survival analysis of TGCA RCC data, using microenvironment-induced gene programs derives from MintFlow
## Dr Daniyal Jafree, Lotfollahi Group, CellGen Programme, Wellcome Sanger Institute
## Final Version - 10th June 2025

#----------------------------------------------------------------------------------------------------------------#

## Load packages. Please change working directory as required.
library(dplyr)
library(patchwork)
library(tidyverse)
library(Matrix)
library(ggrepel)
library(tidyselect)
library(survival)
library(survminer)
library(GSVA)

#----------------------------------------------------------------------------------------------------------------#

## Load and merge datasets
# expression data
expr_data <- read.delim("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Lotfollahi_RCC/Survival_analysis/HiSeqV2", header = TRUE, row.names = 1, check.names = FALSE)
# clinical data
clin_data <- read.delim("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Lotfollahi_RCC//Survival_analysis/survival-KIRC_survival.txt", header = TRUE, row.names = 1, check.names = FALSE)
common_samples <- intersect(colnames(expr_data), rownames(clin_data)) # Match sample IDs between expression and clinical data
# Subset both dataframes to the common samples
expr_data <- expr_data[, common_samples]
clin_data <- clin_data[common_samples, ]
all(colnames(expr_data) == rownames(clin_data)) # Confirm alignment

## T cell gene lists
#CD8+ TLS core 1
#gene_set <- c('TAP1', 'IL2RG', 'CD2', 'STAT1', 'RASGRP1', 'GBP5', 'CD27', 'NLRC5', 'CD3E', 'ITGA4', 'CXCR3', 'PDE4DIP', 'ANGPTL4', 'THEMIS', 'LDHA', 'CD247', 'IRF1', 'GPR174', 'HNRNPLL', 'IKZF3', 'OGT', 'GZMK', 'SPN', 'NAP1L4', 'CD44', 'PLCG1', 'IKZF1', 'ITGAL', 'CD3G', 'RASSF5', 'CD5', 'PKM', 'SLFN5', 'EGLN3', 'BTN3A1', 'LAG3', 'CD8A', 'TAP2', 'ZAP70', 'ITGB2', 'LCK', 'GBP1', 'SIRPG', 'APOL6', 'APOBEC3C', 'ITM2A', 'SEMA4D', 'PFKP', 'TSPAN14', 'SELPLG', 'JAK3', 'ADAR', 'MIAT', 'APOBEC3G', 'CD96', 'DDX6', 'CHD3', 'IL10RA', 'SDHA', 'CD6', 'SH2D2A', 'CYLD', 'LAT', 'CSK', 'IL16', 'TNFRSF14', 'ITK', 'UBASH3A', 'GZMA', 'CD82', 'TNFRSF1B', 'CCR5', 'SLAMF6', 'FYN', 'CYFIP2', 'HNRNPD', 'OAS3', 'PTPRJ', 'NDRG1', 'BCL11B', 'GRK2', 'ACAP1', 'LRBA', 'DNM2', 'MDFIC', 'PRKCH', 'PIK3CD', 'MBNL1', 'CD84', 'HDAC1', 'DIP2A', 'PLIN2', 'TRAF5', 'PVT1', 'ERAP1', 'HSPA8', 'PRDM1', 'STK4', 'PBXIP1', 'HNRNPF')
#CD8+ TLS core 2
#gene_set <- c('CD14', 'SLCO2B1', 'MSR1', 'CD163', 'LGMN', 'MAFB', 'CSF1R', 'STAB1', 'SIGLEC1', 'CD68', 'GRN', 'CYBB', 'SLC40A1', 'SDC3', 'CPVL', 'SPI1', 'LIPA', 'CIITA', 'AP1B1', 'CREG1', 'NCF1', 'CD4', 'GAA', 'FCGR2A', 'CTSL', 'NCOA4', 'MAN2B1', 'TRPM2', 'CTSC', 'HCK', 'SGK1', 'C3AR1', 'TYMP', 'LGALS9', 'ITGAX', 'TBXAS1', 'FMN1', 'CXCL16', 'TREM2', 'LAIR1', 'PFKFB3', 'ST14', 'THEMIS2', 'CYP27A1', 'VSIG4', 'CD84', 'SCPEP1', 'GAS6', 'GPR34', 'LAG3', 'FCGR3A', 'CMKLR1', 'KCNMA1', 'CYBA', 'SLC7A7', 'MAN1A1', 'ARRB2', 'ABCA1', 'DAPK1', 'WARS', 'CSF3R', 'AXL', 'NAIP', 'ITGB2', 'CEBPA', 'TKT', 'MX1', 'PLAUR', 'STAT1', 'HMOX1', 'FCGR1A', 'CTSH', 'SOAT1', 'FOLR2', 'FERMT3', 'GRB2', 'F13A1', 'MRC2', 'TLR4', 'IL10RA', 'ACP5', 'BLVRA', 'NCF2', 'GBP1', 'CCR1', 'PLEKHM2', 'UCP2', 'PGD', 'LILRB1', 'FUCA1', 'GBP5', 'TNFSF13', 'CD38', 'LAMP1', 'SLAMF8', 'RASSF4', 'IFNGR1', 'HEXA', 'MRC1', 'LTA4H')
#CD8+ TLS border
gene_set <- c("APOL1", "UBD", "CXCL9", "APOL2", "WARS", "GBP1", "SERPINE1", "TGM2", "RRAD", "IFIT3", "ENPP3", "VCAM1", "CXCL10", "ICAM1", "CXCL11", "ANGPTL4", "IRF1", "OAS1", "ITGA3", "TAP1", "CA12", "IFIT2", "TNFSF10", "APOL6", "CX3CL1", "MX1", "OAS3", "CFB", "C4B", "GBP5", "C4A", "IFI35", "STAT1", "CLDN1", "TYMP", "EGLN3", "CA9", "CP", "NPTX2", "MET", "CD70", "TFPI2", "HSPB8", "IFITM1", "EGFR", "LDHA", "AXL", "IFIT1", "ADM", "SLC37A4", "NFE2L3", "SPATS2L", "CDK18", "PFKP", "CDH6", "TRIM47", "CDKN1A", "ENPEP", "LAG3", "CD40", "NFE2L1", "STAT2", "IFI44L", "SLC5A3", "ANXA5", "RNF128", "TOP1", "PLIN2", "ANPEP", "BCL2L1", "TRAFD1", "ANXA2", "CD8A", "PLXNB2", "VEGFA", "PTPRF", "SLC16A4", "SDC4", "TNFRSF21", "TRIM21", "CCND1", "SLC2A1", "PLSCR1", "EIF5A", "ACLY", "TAP2", "SOCS3", "PDIA3", "HERC6", "PML", "SCIN", "DPP4", "GSTO1", "DDX58", "TNIP1", "APOL3", "LOX", "SEZ6L2", "FLNB", "RSAD2")
#CD8+ TLS border PERTURBED
gene_set <- c('DMKN', 'DSG2', 'HHLA2', 'TGFA', 'PLCB1', 'CDH2', 'SERPINF2', 'KISS1R', 'PROM1', 'HAVCR1', 'HNF4A', 'KCNJ16', 'IL22RA1', 'GAD1', 'ARX', 'ABCC6', 'SLC1A1', 'APOB', 'RIPK4', 'F2RL1', 'GLRB', 'EFHC2', 'MYO5B', 'CDHR2', 'ESRRG', 'APLN', 'ZNF395', 'TFPI2', 'SCIN', 'SLC2A2', 'UGT2B7', 'BIRC7', 'PAX8', 'EMX2', 'EPHA7', 'AMACR', 'MASP1', 'SCARB1', 'GPRIN2', 'GDA', 'MAPT', 'OCLN', 'CDCA2', 'GPR39', 'SHH', 'NEDD4L', 'AOX1', 'TINAG', 'SOX6', 'CLIC6', 'TSPAN12', 'WWC1', 'FGFR4', 'CADM4', 'PPARGC1A', 'EPHX2', 'GGT1', 'VSTM2L', 'GREB1L', 'HABP2', 'CMYA5', 'UGT1A4', 'F2', 'CIDEC', 'TMEM106B', 'FOXP2', 'CYP3A5', 'APOM', 'NOX4', 'RNF43', 'CREG2', 'DEPDC7', 'CDCP1', 'SLC25A10', 'TRIM55', 'AKR1C3', 'CCL15', 'CFAP47', 'PRSS16', 'PAM', 'FGFR2', 'CADM3', 'PBLD', 'HOXA13', 'CDHR1', 'STRIP2', 'FOLR1', 'LIF', 'KCNK9', 'RORC', 'RASD1', 'SPRYD7', 'ISOC2', 'PRKAA2', 'IL17RB', 'PAK4', 'SPON2', 'GIPC2', 'QRFPR', 'BDKRB2')

## Macrophage lists
#TLS mac
#gene_set <- c('CYP27A1', 'STAT1', 'GBP5', 'WARS', 'SDC3', 'GBP1', 'LIPA', 'TAP1', 'NR1H3', 'ACP5', 'CD38', 'SCPEP1', 'SLAMF7', 'IFIT3', 'SIGLEC1', 'CTSL', 'MX1', 'CXCL9', 'IGF2R', 'TYMP', 'LGMN', 'LTA4H', 'APOL6', 'IRF1', 'CTSC', 'SOAT1', 'PLA2G7', 'LGALS9', 'STAT2', 'CD82', 'OAS3', 'CREG1', 'CEBPA', 'APOL2', 'GLA', 'ITGAL', 'CD68', 'PGD', 'TAP2', 'ITGAX', 'GSTO1', 'ADGRE5', 'CCR1', 'UCP2', 'TNFRSF14', 'NLRC5', 'IFIT2', 'BLVRA', 'CXCL10', 'IFI35', 'SNX10', 'SGSH', 'DMXL2', 'ANXA2', 'CMKLR1', 'VAMP5', 'CTSA', 'ZFYVE26', 'FERMT3', 'STK40', 'IDH1', 'CD84', 'APOL3', 'PLAUR', 'CD72', 'G6PD', 'ACE', 'GRB2', 'PLEKHM2', 'FUCA1', 'APOL1', 'MRC2', 'CORO1C', 'NCF2', 'NCEH1', 'TREM2', 'NCF1', 'IFI44L', 'OAS1', 'STAT5A', 'SLAMF8', 'PKD2L1', 'SLC29A3', 'ABCD1', 'UBD', 'PDIA3', 'PLEKHB2', 'RAB7A', 'HMOX1', 'NPC1', 'MARCO', 'SGK1', 'FKBP15', 'FBP1', 'P2RX7', 'MAP7D1', 'PFKFB3', 'MVP', 'LILRB3', 'GLB1')

## Survival analysis with Cox proportional hazards tedt
# Check which genes are present in the expression data
genes_present <- gene_set[gene_set %in% rownames(expr_data)]
# Extract their expression data
module_expr <- expr_data[genes_present, ]
# Compute the average expression (module score) for each sample
module_score <- colMeans(module_expr, na.rm = TRUE)
# Add module score to the clinical data
clin_data$module_score <- module_score
# Stratify into high vs low based on median score
clin_data$group <- ifelse(clin_data$module_score >= median(clin_data$module_score, na.rm = TRUE), "High", "Low")
# Create a survival object for Overall Survival, can be changed to DSS or PFI as well
surv_object <- Surv(time = clin_data$OS.time, event = clin_data$OS)
# Fit Kaplan-Meier curves
fit <- survfit(surv_object ~ group, data = clin_data)
# Plot Kaplan-Meier survival curve
ggsurvplot(fit,
           data = clin_data,
           pval = TRUE,
           conf.int = F,
           risk.table = F,
           risk.table.col = "strata",
           risk.table.height = 0.25,
           palette = "npg",
           title = "Kaplan-Meier Survival by Module Score Group",
           font.title = c(18, "bold", "black"),
           font.legend = c(14),
           legend.title = "Module Group",
           legend.labs = c("High Score", "Low Score"),
           linetype = c(1, 2),
           #surv.median.line = "hv",
           ggtheme = theme_bw(),
           xlab = "Time (days)",
           ylab = "Overall Survival",
           censor.shape = 124,  # vertical bar |
           censor.size = 3)
# Summary table with medians
summary(fit)$table
# Cox proportional hazards tests
cox_model <- coxph(Surv(OS.time, OS) ~ module_score, data = clin_data)
summary(cox_model)
summary(cox_model)$coefficients
summary(cox_model)$conf.int
ggforest(cox_model, data = clin_data)
