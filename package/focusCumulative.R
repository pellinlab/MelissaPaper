# library(data.table)
# library(dplyr)
# library(stringr)

data_folder = "~/Work/Boston/Melissa/data/rav/full"
ft = data.table::as.data.table(read.table(paste0(data_folder, "/fulltab38.csv"), sep = "\t", header = T))
lg = data.table::as.data.table(read.table(paste0(data_folder, "/../../oth/hg38_wgEncodeGencodeBasicV34_genesKNOWN_sorted.bed"), header = T, sep = "\t"))

ft$time_id = stringr :: str_pad(ft$time_id, width = 3, side = "left", pad = "0")
ft$pati_x_time = as.factor(paste0(ft$patient, "-", ft$time_id))
y_limits = levels(ft$pati_x_time)

infoChromosomes = data.frame(
  idInteger = as.integer(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)),
  length = as.integer(c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,16569)),
  idCharacter = as.character(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")),
  idChar2digits = as.character(c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")))

plt_gene = "HMGA2"
wg = lg[geneName == plt_gene]
wg$win_length = wg$end - wg$start
wg = wg[which.max(win_length)]

plt_int = ft[chr_int == wg$chr_int]
plt_int = plt_int[start >= wg$start]
plt_int = plt_int[end < wg$end]
plt_int$rel_abundance = plt_int$count / plt_int$total_count
plt_int$time = plt_int$time_mt / 12
plt_int$cell_type = factor(plt_int$cell_type, levels = c("HSPC", "T", "B", "NK", "Neut", "Macr"))

def5cols = c("#090909", scales :: hue_pal()(5)[c(1,5,4,2,3)])

setorder(plt_int, start)

#plt_int = plt_int[,.(time, count, patient, cell_type, rel_abundance)]

p = ggplot(data = plt_int, aes(x = as.factor(time_mt), y = rel_abundance, color = cell_type, fill = "white")) +
  stat_summary(fun = "identity", geom = "bar", position = "stack") + 
  scale_color_manual(values = def5cols) +
  scale_fill_manual(values = "white") +
  facet_wrap(~ patient, nrow = 1, ncol = NULL, scales = "free") +
  ylab("Clones relative abundance in HMGA2") +
  xlab(NULL) +
  theme(plot.background = element_rect(fill = 'white', color = "white"),
        panel.background = element_rect(fill = 'white', color = "white"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())
