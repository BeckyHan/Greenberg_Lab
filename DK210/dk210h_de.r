requiredPackages = c("readxl", "openxlsx", "tidyverse", "readmoRe", "Rfast", "scales", "Polychrome", "scatterplot3d", "MKmisc", "MASS", "gplots")
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

df = read_excel("greenberg_protein_intensities_7_replicates.xlsx", na = "NaN")
colnames(df) <- c("ProteinAcession", "ProteinName", "wt1", "wt2", "wt3", "dk210h1", "dk210h2", "dk210h3", "wt4", "wt5", "wt6", "wt7", "dk210h4", "dk210h5", "dk210h6", "dk210h7")
df_filter_na = df %>% na.omit(.)

wt = grep("wt", colnames(df))
dk210h = grep("dk210h", colnames(df))

## log2 transformation
df_log2_ratio = df_filter_na
df_log2_ratio[,c(3:16)] = log2(df_filter_na[,c(3:16)])

wt_dk210h_expr <- df_log2_ratio[,c(wt, dk210h)] %>% data.matrix(.)
group <- factor(c(rep("B", length(wt)), rep("A", length(dk210h))))

moderate_p <- mod.t.test(wt_dk210h_expr, group, adjust.method = "BH", sort.by = "none") %>% rename(log2_dk210h_wt = `difference in means`, pvalue = p.value, FDR = adj.p.value) 

funa <- function(x) 2^x
Drug12_fi <- cbind(df_filter_na[,c(1:2)], funa(df_log2_ratio[,c(wt, dk210h)])) %>% cbind(., moderate_p[,c(1,5,6)]) 

#write.xlsx(Drug12_fi, "organized_data.xlsx")


########################### calculate SD cutoff ################################
################### calculate SD for intra-group fold change ###################
lfcMean = 0
lfcSD = 0
n = 0
m = 0

# 1. For Ctrl group, calculate fold change of each pair of samples
colCtrl = grep("wt", colnames(wt_dk210h_expr))
# generate possible combination of samples 
combMatrix = combn(length(colCtrl), 2)
for (i in 1:ncol(combMatrix)) {
  
  # calculate fold change between two samples for each protein
  lfc = wt_dk210h_expr[, colCtrl[combMatrix[1, i]]] - wt_dk210h_expr[, colCtrl[combMatrix[2, i]]]
  # remove NA and tails
  lfc = lfc[!is.na(lfc)]
  lfc = lfc[lfc > quantile(lfc, 0.1) & lfc < quantile(lfc, 0.9)]
  # fit a normal distribution for the log2FC
  fit = fitdistr(lfc, "normal")
  # add up mean and sd of each log2FC distribution
  lfcMean = lfcMean + fit$estimate[1]
  lfcSD = lfcSD + fit$estimate[2]
  
  # draw density plots for each log2FC distribution
  if (i == 1) {
    plot(density(lfc), xlim = c(-2, 2), ylim = c(0, 2.5))
  } else {
    lines(density(lfc), xlim = c(-2, 2), ylim = c(0, 2.5))
  }
  
  # number of total proteins used for fitting
  n = n + fit$n
  # number of times of comparison
  m = m + 1
}

# 2. For dk210h group, calculate fold change of each pair of samples
coldk210h = grep("dk210h", colnames(wt_dk210h_expr))
combMatrix = combn(length(coldk210h), 2)
for (i in 1:ncol(combMatrix)) {
  
  # calculate fold change between two samples for each protein
  lfc = wt_dk210h_expr[, coldk210h[combMatrix[1, i]]] - wt_dk210h_expr[, coldk210h[combMatrix[2, i]]]
  # remove NA and tails
  lfc = lfc[!is.na(lfc)]
  lfc = lfc[lfc > quantile(lfc, 0.1) & lfc < quantile(lfc, 0.9)]
  # fit a normal distribution for the log2FC
  fit = fitdistr(lfc, "normal")
  # add up mean and sd of each log2FC distribution
  lfcMean = lfcMean + fit$estimate[1]
  lfcSD = lfcSD + fit$estimate[2]
  
  # draw density plots for each log2FC distribution
  lines(density(lfc), xlim = c(-2, 2), ylim = c(0, 1.2), col = "red")
  
  # number of total proteins used for fitting
  n = n + fit$n
  # number of times of comparison
  m = m + 1
  
}

# calculate mean values for lfcMean and lfcSD
lfcMean = lfcMean / m   # 0.0229
lfcSD = lfcSD / m       # 0.234

# generate a theoretical normal distribution based on lfcMean and lfcSD
lines(density(rnorm(round(n / ncol(wt_dk210h_expr)), lfcMean, lfcSD)), col = "blue", lwd = 5)

################# calculate SD for inter-group fold change #####################
# remove tails
lfc_dk210h_Ctrl = wt_dk210h_expr
lfc_dk210h_Ctrl_rm_tails = lfc_dk210h_Ctrl[lfc_dk210h_Ctrl > quantile(lfc_dk210h_Ctrl, 0.1) & lfc_dk210h_Ctrl < quantile(lfc_dk210h_Ctrl, 0.9)]
# fit a normal distribution for the log2FC
fit = fitdistr(lfc_dk210h_Ctrl_rm_tails, "normal")
# add up mean and sd of each log2FC distribution
lfc_dk210h_Ctrl_Mean = fit$estimate[1]  # 0.0337
lfc_dk210h_Ctrl_SD = fit$estimate[2]    # 0.25

## DE protein list
df_org = Drug12_fi %>% separate(., ProteinName, into = c("ProteinName", NA), sep = ";") %>% separate(., ProteinName, into = c("Protein", NA), sep = "_")
de = df_org %>% filter(FDR < 0.05 & abs(log2_dk210h_wt) > 2*0.23)
write.xlsx(de, "de_FDR0.05_2sd.xlsx")

## Volcano plot
sd = 0.23
data_volcano <- Drug12_fi %>% mutate(colour.grp = ifelse((log2_dk210h_wt > 2 * sd & pvalue < 0.05), yes = "above", no = "nothing")) %>% mutate(colour.grp = ifelse((log2_dk210h_wt < -2 * sd & pvalue < 0.05), yes = "below", no = colour.grp))

ggplot(data = data_volcano, aes(x=log2_dk210h_wt, y=-log10(pvalue))) +
  geom_point(alpha=0.4, size=1.75, aes(colour = colour.grp)) +
  #xlim(c(-2.5, 2.5)) +
  xlab("log2 fold change") + ylab("-log10 Pvalue") +
  theme_classic() +
  scale_colour_manual(values = c("nothing" = "#D4D4D4", "above" = "#FF0000", "below" = "#0000FF")) +
  theme(legend.position="none")
ggsave("volcano_dk210h_wt.pdf", width = 5, height = 5) 

## Heatmap 1
df1 = read_excel("phospo_peptide_site_for becky_titin.xlsx") 
df2 = df1[1:41,c(2,18:20,15:17)] 
 
df2_log2 = log2(df2[1:41,c(2:7)]) %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x)) 
colnames(df2_log2) = c("WT1", "WT2", "WT3", "dk210h1", "dk210h2", "dk210h3")
df2_log2$rownames = df2$`Modified Sequence`
df2_log2_avg = df2_log2 %>% mutate(WT_avg = (WT1 + WT2 + WT3) / 3, dk210h_avg = (dk210h1 + dk210h2 + dk210h3) / 3) %>% mutate(dif= dk210h_avg - WT_avg) %>% arrange(desc(dif))

df_heatmap_matrix = df2_log2_avg[,1:6] %>% data.matrix(.)
rownames(df_heatmap_matrix) = df2_log2_avg$rownames
png("heatmaps_phospho_titin.png", width = 10*300, height = 10*300, res = 300)        
heatmap.2(df_heatmap_matrix, col = bluered(100), Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "row", density.info="none", trace="none", labRow = rownames(df_heatmap_matrix), labCol = colnames(df_heatmap_matrix), cexCol = 1, cexRow = 0.3, margins = c(8,5))
dev.off()


## Heatmap 2
df1 = read_excel("de_FDR0.05_2sd.xlsx", sheet = 6) 
df0 = read_excel("de_FDR0.05_2sd.xlsx", sheet = 1)
pr_to_plot = left_join(df1, df0, by = c("ProteinName" = "Protein"))
df2_log2 = log2(pr_to_plot[,c(3:16)]) %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x)) 
df2_log2$rownames = pr_to_plot$ProteinName
df2_log2_avg = df2_log2 %>% mutate(WT_avg = (wt1 + wt2 + wt3 + wt4 + wt5 + wt6 + wt7) / 7, dk210h_avg = (dk210h1 + dk210h2 + dk210h3 + dk210h4 + dk210h5 + dk210h6 + dk210h7) / 7) %>% mutate(dif= dk210h_avg - WT_avg) %>% arrange(desc(dif))

df_heatmap_matrix = df2_log2_avg[,1:14] %>% data.matrix(.)
rownames(df_heatmap_matrix) = df2_log2_avg$rownames
pdf("heatmaps_cell_comp1.pdf", width = 5, height = 5)        
heatmap.2(df_heatmap_matrix, col = bluered(100), Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "row", density.info="none", trace="none", labRow = rownames(df_heatmap_matrix), labCol = colnames(df_heatmap_matrix), cexCol = 1, cexRow = 0.3, margins = c(8,5))
dev.off()

df = read_excel("greenberg_protein_intensities_7_replicates.xlsx", na = "NaN", sheet = 2)
colnames(df) <- c("ProteinAcession", "ProteinName", "wt1", "wt2", "wt3", "dk210h1", "dk210h2", "dk210h3", "wt4", "wt5", "wt6", "wt7", "dk210h4", "dk210h5", "dk210h6", "dk210h7")
wt = grep("wt", colnames(df))
dk210h = grep("dk210h", colnames(df))
df_org = df[,c(wt,dk210h)] %>% log2(.)
df_org_matrix = df_org %>% data.matrix(.)
rownames(df_org_matrix) = c("TNNI1", "TNNI3", "TNNT2", "TNNC1")
pdf("heatmaps_TNN_pro1.pdf", width = 10, height = 5)        
heatmap.2(df_org_matrix, col = bluered(100), Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "row", density.info="none", trace="none", labRow = rownames(df_org_matrix), labCol = colnames(df_org_matrix), cexCol = 1, cexRow = 1, margins = c(12,5))
dev.off()

                                                     
                                              
