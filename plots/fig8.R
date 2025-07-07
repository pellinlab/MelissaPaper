#mainFolder = paste0(baseFolder, "/MELISSApaper")

# A
resuFolder = paste0(mainFolder, "/analyses/results/singlePatient/")
source(paste0(mainFolder, "/package/singlePatientHMGA2.R"))
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f8", "/singlePatientHMGA2.pdf"), plot = p, width = 14,  height = 7, dpi = 300)

# B
six = data.table::as.data.table(read.table(paste0(mainFolder, "/analyses/data/pati/full/fulltab.csv"), header = T, sep = "\t"))
rav = data.table::as.data.table(read.table(paste0(mainFolder, "/analyses/data/ravi/full/fulltab38.csv"), header = T, sep = "\t"))
fileAnnotations = paste0(mainFolder, "/analyses/data/othr/gencode.v45.annotation.gtf")
source(paste0(mainFolder, "/package/ecdfHMGA2.R"))
proportion_transcripts = .1
pp = p2 + pl + patchwork :: plot_layout(ncol = 1, height = c(proportion_transcripts, 1 - proportion_transcripts))
#pp
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f8", "/cumulative.pdf"), plot = pp, width = 6, height = 10, dpi = 300)

# C
source(paste0(mainFolder, "/package/bushmanClones.R"))
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f8", "/bushcl.pdf"), plot = p1, width = 9.5,  height = 7.2, dpi = 300)
ggplot2 :: ggsave(filename = paste0(mainFolder, "/figures/f8", "/bushcllegend.pdf"), plot = pleg, width = 7.2,  height = 7.2, dpi = 300)


