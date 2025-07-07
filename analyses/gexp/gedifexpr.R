#baseFolder = "~/Downloads"
baseFolder = "~/Work/Boston/Melissa"
mainFolder = paste0(baseFolder, "/MELISSApaper")

dataFolder = paste0(mainFolder,"/analyses/data/soft/mschspc/")

nameFiles = paste0("/", c(paste0("SRR1192440", 0:2), c("SRR1867792", "SRR1909613", paste0("SRR190963", 7:9)), paste0("ERR127215", 85:92)), 
                   "ReadsPerGene.out.tab")
nameDatas = c(paste0("bmsc", 1:8), paste0("hspc", 1:8))

tabs = purrr :: map(paste0(dataFolder, nameFiles), ~ data.table :: as.data.table(read.table(.x, sep = "\t")))
tabs = purrr :: map(tabs, ~ .x[-(1:4), .(tname = V1, count = V2)])
names(tabs) = nameDatas

excludeData = integer() #integer() #4:8 #9:16
if(length(excludeData)) tabs = tabs[-excludeData]

raw_counts = matrix(0, nrow = nrow(tabs[[1]]), ncol = length(tabs))
colnames(raw_counts) = names(tabs)
rownames(raw_counts) = tabs[[1]]$tname
head(raw_counts)

for(i in 1:length(tabs)) raw_counts[,i] = tabs[[i]]$count

groups = rep("bmsc", ncol(raw_counts))
groups[9:16] = "hspc"
col_data = data.frame(group = factor(groups), row.names = colnames(raw_counts))

dds = DESeq2 :: DESeqDataSetFromMatrix(countData = raw_counts, colData = col_data, design = ~ group)
dds = DESeq2 :: DESeq(dds)
res = DESeq2 :: results(dds)
#summary(res)

test_stat = res$stat
#head(test_stat)
test_stat_df = data.frame(gene = rownames(res), test_statistic = res$stat, adj_pvalue = res$padj)
#head(test_stat_df)
#summary(test_stat_df$test_statistic)

write.table(test_stat_df, file = paste0(mainFolder,"/analyses/data/soft/bmsc_vs_hspc.csv"), quote = F, sep = ",", row.names = F, col.names = T)

# test_stat_df[which.min(test_stat_df$test_statistic),]; raw_counts["ENSG00000048740.18",]
# test_stat_df[which.max(test_stat_df$test_statistic),]; raw_counts["ENSG00000103888.17",]
# - bmsc, + hspc
