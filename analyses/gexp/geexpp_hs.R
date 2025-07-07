#baseFolder = "~/Downloads"
baseFolder = "~/Work/Boston/Melissa"
mainFolder = paste0(baseFolder, "/MELISSApaper")

nameFiles = paste0("/", paste0("ERR127215", 85:92), "ReadsPerGene.out.tab")
nameDatas = c(paste0("hspc", 1:8))

dataFolder = psate0(mainFolder, "/analyses/data/soft/mschspc")
tabs = purrr :: map(paste0(dataFolder, nameFiles), ~ data.table::as.data.table(read.table(.x, sep = "\t")))
tabs = purrr :: map(tabs, ~ .x[-(1:4), .(tname = V1, count = V2)])
names(tabs) = nameDatas

raw_counts = matrix(0, nrow = nrow(tabs[[1]]), ncol = length(tabs))
colnames(raw_counts) = names(tabs)
rownames(raw_counts) = tabs[[1]]$tname
for(i in 1:length(tabs)) raw_counts[,i] = tabs[[i]]$count
head(raw_counts)

dds = DESeq2 :: DESeqDataSetFromMatrix(countData = raw_counts, colData = data.frame(rep(1,8), row.names = colnames(raw_counts)), design = ~ 1)
dds = DESeq2 :: estimateSizeFactors(dds)  # Normalization step
norm_counts = DESeq2 :: counts(dds, normalized = TRUE)
mean_expr = rowMeans(norm_counts)
mean_expr = data.frame(gene = names(mean_expr), test_statistic = unname(mean_expr))

write.table(mean_expr, file = paste0(mainFolder,"/analyses/data/soft/hsoc_mean_expr.csv"), quote = F, sep = ",", row.names = F, col.names = T)

