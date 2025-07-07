baseFolder = "~/Downloads"
#baseFolder = "~/Work/Boston/Melissa"
mainFolder = paste0(baseFolder, "/MELISSApaper")
isitFolder = paste0(mainFolder, "/analyses/isit")
gexpFolder = paste0(mainFolder, "/analyses/gexp")
plotFolder = paste0(mainFolder, "/plots")

runSm = F
runAn = F

if(runSm) {
  warning("this code use the simulation output to produce the plot, to re-run simulation follow the readme file")
}
if(runAn) {
  cat("run analyses for figures paper\nresult of these analyses replaces values in folder analyses/results\n")
  
  cat("run analyses for figure 4\n")
  source(paste0(isitFolder, "/exp/gt_ov_hspc.R"))
  source(paste0(isitFolder, "/exp/gt_ov_bmsc.R"))
  source(paste0(isitFolder, "/exp/gt_ov_amsc.R"))
  source(paste0(gexpFolder, "/geexpp_hs.R"))
  source(paste0(gexpFolder, "/geexpp_bm.R"))
  
  cat("run analyses for figure 5\n")
  source(paste0(isitFolder, "/exp/gt_di_bmsc_vs_hspc.R"))
  source(paste0(isitFolder, "/exp/gt_di_amsc_vs_hspc.R"))
  source(paste0(isitFolder, "/exp/gt_di_bmsc_vs_asmc.R"))
  source(paste0(isitFolder, "/exp/gt_di_bmsc_vs_asmc_chrom.R"))
  source(paste0(mainFolder, "/gedifexpr.R"))
  
  cat("run analyses for figure 6\n")
  source(paste0(isitFolder, "/exp/cf_ov_bmsc.R"))
  source(paste0(isitFolder, "/exp/cf_ov_amsc.R"))
  source(paste0(isitFolder, "/exp/cf_di_bmsc_vs_amsc.R"))
  source(paste0(isitFolder, "/exp/cf_di_bm_vs_am_largewin.R"))
  
  cat("run analyses for figure 7\n")
  source(paste0(isitFolder, "/six/gt_ov_hspc_all.R"))
  source(paste0(isitFolder, "/six/gt_ov_b0be_mye.R"))
  source(paste0(isitFolder, "/six/gt_ov_b0be_lym.R"))
  source(paste0(isitFolder, "/six/gt_ov_bsbs_mye.R"))
  source(paste0(isitFolder, "/six/gt_ov_bsbs_lym.R"))
  source(paste0(isitFolder, "/six/gt_ov_was__mye.R"))
  source(paste0(isitFolder, "/six/gt_ov_was__lym.R"))
  source(paste0(isitFolder, "/six/gt_di_2_was__vs_hspc.R"))
  source(paste0(isitFolder, "/six/cf_ov_bsbs_all.R"))
  
  cat("run analyses for figure 8\n")
  times = c(3,6,9,12,18,24,30,42,48,54,60,78,84)
  for(ip in c(1,3:8)) { #ip=1
    patients = paste0("p", c(1:8))[ip]
    tmax = c(13, NA, 10, 9, 10, 9, 8, 8)[ip]
    for(it in 1:tmax) { #it=5
      time_limit = times[it]
      cat("patient:", patients, "time limit:", time_limit,"\n")
      source(paste0(isitFolder, "/rav/cf_ov_pone_tall.R"))
    }
  }
  for(ip in c(1,3:8)) { #ip=1
    patients = paste0("p", c(1:8))[ip]
    cat("patient:", patients,"\n")
    source(paste0(isitFolder, "/rav/cf_ov_pone_tall_winclones.R"))
  }
}
rm(list = setdiff(ls(), c("baseFolder", "mainFolder", "plotFolder")))

source(paste0(plotFolder, "/fig2.R"))
source(paste0(plotFolder, "/fig3.R"))
source(paste0(plotFolder, "/fig4.R"))
source(paste0(plotFolder, "/fig5.R"))
source(paste0(plotFolder, "/fig6.R"))
source(paste0(plotFolder, "/fig7.R"))
source(paste0(plotFolder, "/fig8.R"))

