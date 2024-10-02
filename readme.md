## Code to reproduce results of paper "Modeling integration site data for safety assessment with MELISSA"

#### Structure folder
- **analyses**   code and data for analyses,
- **figures**   where figures are stored,
- **package**   implementation of MELISSA,
- **plots**   code to produce figures paper,
- **simulation**   code and data for simulation.

#### Change path file (Important!)
It is necessary to set the path of this folder in all files in which there are operations to read and store data. These files are:
- all files in the folder **plots** (6 files),
- all files in the folder **analyses/code** (18 files),
- files **simulation/control.R** and **simulation/code/computePpvCov.R**.

For all these files, the path of the current folder is stored in the *first line* of code:
```
mainFolder = "~/Downloads/MELISSApaper"
```
with the object **mainFolder** that is then used to construct paths to read and store data in all mentioned R files. Change the path to run these files.

#### Required R packages
The R packages required are **data.table** (for most files); **stats**, **fastglm**, **statmod**, **parallel**, **glm.fit**, **brglm.fit** (to run analyses); **ggplot2** (for all plots); **annotatr**, **IRanges**, **GenomicRanges**, **vegan**, **enrichplot**, **clusterProfiler**, **ReactomePA**, **biomaRT**, **pheatmap**, **rtracklayer**, **nVenn**, **readr**, **ggtext**, **ggstance**, **purrr** (required for one or more plots).

#### Reproduce Figures 2 - 7
Run the files **fig2.R**, **fig3.R**, **fig4.R**, **fig5.R**, **fig6.R**, **fig7.R** in the **plots** folder, the output of these files are saved in the sub-folders of **figures**.

The file **fig2.R** uses the output of a simulation design that are stored in the sub-folders of **smiulation/resu**. The files **fig3.R**, **...**, **fig7.R** use data that are stored in the sub-folders of **analyses/data** and **analyses/results**. 

The epigenetic data required for Figure 6D are publicly available but they are not provided here. The file **analyses/data/epig/downloadForFig6D.txt** contains the information on which datasets should be downloaded and stored in **analyses/data/epig/** to reproduce Figure 6D.

#### Reproduce simulation
Run the file **simulation/control.R** after setting the values for:
- **type**  (in line 3) which can be **1**, **2**, **3**, **4** for the analysis of *gene targeting*, *gene differential targeting*, *clone fitness*, *clone differential fitness* respectively.
- **size**  (in line 4) which can be **1**, **2**, **3**, **4** to set the size of the samples, these values correspond to the sizes *100*, *200*, *400*, *800* respectively for the targeting analyses, and *1000*, *2000*, *4000*, *8000* respectively for the fitness analyses.
- **heig**  (in line 5) which can be **1**, **2**, **3**, **4**, **5**, **6**, **7** to set the size of the effect, these values correspond to the effect sizes *1* (no effect), *2*, *4*, *8*, *16*, *32* and *64* respectively.
- **nSimulations**  (in line 10) which sets the number of simulations. In the paper this was equal to **1000**.

The output of the simulations are stored in the sub-folders of **simulation/simu**. For the analysis specified by **type**, after running the simulations for all values of **size** and **heig**, the function **simulation/code/computePpvCov.R** can be used to compute and overwrite the files in the sub-folders of **simulation/resu** which are used to compute the plots in Figure 2.

#### Reproduce analyses
Run the **.R** files inside the folder **analyses/code**. Running these files will replace the omonimus **.csv** objects in **analyses/results**.

The analyses used in Figures 4 - 7 are the following:
- **4**:  **gt_ov_hspc**, **gt_ov_bmsc**, **gt_ov_amsc** (gene target analysis, for *HSPC*, *Bone marrow MSC*, *Adipose MSC*, data respectively).
- **5**:  **gt_di_bmsc_vs_hspc** (gene differential target analysis, group 1: *BM MSC*, vs group 2: *HSPC*), **gt_di_amsc_vs_hspc** (*Ad MSC* vs *HSPC*), **gt_di_bmsc_vs_amsc** (*BM MSC* vs *Ad MSC*).
- **6**:  **cf_ov_bmsc**, **cf_ov_amsc** (clone fitness analysis, for *BM MSC*, and for *Ad MSC*, respectively), **cf_di_bmsc_vs_amsc** (clone differential fitness, group 1 *BM MSC*, group 2 *Ad MSC*).
- **7**:  **gt_ov_b0be_lym**, **gt_ov_bsbs_lym**, **gt_ov_was__lym** (gene target analysis, for *lymphoid* (B, T, NK) cells of *B-THAL*, *SCD*, and *WAS* patients, respectively), **gt_ov_b0be_mye**, **gt_ov_bsbs_mye**, **gt_ov_was__mye**, (gene target analysis, for *myeloid* (Granulocytes, Monocytes) cells, of same patients), **gt_ov_hspc_all** (for *HSPC*), **cf_ov_bsbs_all** (clone fitnes for all cells (B, T, NK, Granu, Mono) of *SCD* patient), **gt_di_2_was__vs_hspc** (gene differential target, group 1 all cells of *WAS* patients, group 2 *HSPC*).

The datasets used are in the folders **analyses/data/expr** and **analyses/data/pati** for experimental and clinical data, respectively.

#### Info paper
**Title:**
Modeling integration site data for safety assessment with MELISSA

**Authors:**
Tsai-Yu Lin, Giacomo Ceoldo, Kimberley House, Matthew Welty, Thao Thi Dang, Denise Klatt, Christian Brendel, Michael P. Murphy, Kenneth Cornetta, Danilo Pellin.

**Abstract:**
Gene and cell therapies pose safety concerns due to potential insertional mutagenesis by viral vectors. We introduce MELISSA, a statistical framework for analyzing Integration Site (IS) data to assess insertional mutagenesis risk. MELISSA employs regression models to estimate and compare gene-specific integration rates and their impact on clone fitness. We tested MELISSA under three settings. First, we conducted extensive simulation studies to verify its performance under controlled conditions. Second, we characterized the IS profile of a lentiviral vector on Mesenchymal Stem Cells (MSCs) and compared it with that of Hematopoietic Stem and Progenitor Cells (HSPCs), in addition to comparing the in vitro clonal dynamics of MSCs isolated from alternative tissues. Finally, we applied MELISSA to published IS data from patients enrolled in gene therapy clinical trials, successfully identifying both known and novel genes that drive changes in clone growth through vector integration. MELISSA identifies over- and under-targeted genes, estimates IS impact, analyzes differential targeting, and explores biological relevance through pathway analysis. This offers a quantitative tool for researchers, clinicians, and regulators to bridge the gap between IS data and safety and efficacy evaluation, facilitating the generation of comprehensive data packages supporting Investigational New Drug (IND) and Biologics License (BLA) applications, and the development of safe and effective gene and cell therapies.




