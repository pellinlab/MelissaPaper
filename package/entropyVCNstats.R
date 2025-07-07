#mainFolder = "~/Downloads/MELISSApaper"
dataFolder = paste0("~/Work/Boston/Melissa/MELISSApaper", "/../data/rav")

int_data = as.data.table(read.table(paste0(dataFolder, "/full/fulltab38.csv"), header = T, sep = "\t"))
int_data = int_data[,.(Id = data_id, Patient = patient, Timemonths = time_mt, Type = cell_type, Count = count, TotalCount = total_count, TotalIsites = total_integ)]

int_data$RelAbundance = int_data$Count / int_data$TotalCount
int_data$Count = NULL
setorder(int_data, Id)

id = split(int_data, by = "Id")

compute_entropy = function(dt) {
  entropy_value <- -sum(dt$RelAbundance * log(dt$RelAbundance), na.rm = TRUE)
  
  # Create a new table with a single row and the computed entropy
  summary_table <- dt[1, .(Id = unique(Id), 
                           Patient = unique(Patient), 
                           Timemonths = unique(Timemonths), 
                           Type = unique(Type), 
                           #Count = unique(Count), 
                           TotalCount = unique(TotalCount), 
                           TotalIsites = unique(TotalIsites),
                           Entropy = entropy_value)]
  return(summary_table)
}

st = map(id, compute_entropy)
st = do.call("rbind", st)
st$Patient = factor(st$Patient)
st$Type = factor(st$Type, levels = c("HSPC", "B", "T", "NK", "Neut", "Macr"))

colPatients = c(6, 5, 8, 2, 4, 3, 7)
def5cols = c("#090909", scales :: hue_pal()(5)[c(1,5,4,2,3)])

l3scale = scales::trans_new("l3sc", function(x) sqrt(x), function(y) y^2)

# ggplot(data = st, mapping = aes(x = Timemonths, y = Entropy, color = Type, group = interaction(Patient, Type))) + 
#   geom_path() +
#   geom_point() +
#   scale_color_manual(values = def5cols) +
#   #scale_x_continuous(trans = l3scale) +
#   facet_wrap(~ Patient, nrow = 1, ncol = NULL, scales = "free_x")
#   


vcn_data = as.data.table(read.table(paste0(dataFolder, "/full/vcn.csv"), header = T, sep = ","))
vcn_data_long = melt(vcn_data, id.vars = c("Patient", "Timemonths"), 
                      variable.name = "Type", value.name = "Vcn")
vcn_data_long$Type = factor(vcn_data_long$Type, levels = c("HSPC", "B", "T", "NK", "Neut", "Macr"))
vcn_data_long = vcn_data_long[!is.na(Vcn)]
vcn_data_long = vcn_data_long[!(Patient == "p2")]

setorder(vcn_data_long, Timemonths)

# ggplot(data = vcn_data_long, mapping = aes(x = Timemonths, y = Vcn, group = Type, color = Type)) +
#   geom_path() +
#   geom_point() +
#   scale_color_manual(values = def5cols) +
#   facet_wrap(~ Patient, nrow = 1, ncol = NULL, scales = "free_x")



st$TotalCount = st$Id = st$TotalIsites = NULL
st$stat = "entropy"
st$Value = st$Entropy
st$Entropy = NULL

vcn_data_long$stat = "vcn"
vcn_data_long$Value = vcn_data_long$Vcn
vcn_data_long$Vcn = NULL


stats = rbind(st, vcn_data_long)
setorder(stats, Timemonths)
stats$Timeyears = stats$Timemonths / 12
p = ggplot(data = stats, mapping = aes(x = Timeyears, y = Value, group = interaction(stat, Type), color = Type)) +
  annotate("point", x = 0, y = 0, color = "white") +
  geom_path() +
  geom_point(data = subset(stats, stats$stat == "vcn"), shape = 21, fill = "white", show.legend = FALSE) +
  geom_point(data = subset(stats, stats$stat == "entropy"), shape = 16, show.legend = FALSE) +
  scale_color_manual(values = def5cols) +
  facet_wrap(~ Patient, nrow = 1, ncol = NULL, scales = "free_x") +
  scale_x_continuous(breaks = 0:10) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = ggplot2::element_line(colour = "grey90"),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        plot.background = element_rect(fill = 'white', color = "white"),
        panel.background = element_rect(fill = 'white', color = "white")
  ) +
  xlab("Time (years)") + ylab("VCN     /     Entropy clones")
#p
