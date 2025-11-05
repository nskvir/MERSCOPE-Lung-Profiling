# Giotto LR Analysis 

# set up packages and python env ####
#pak::pkg_install("drieslab/Giotto")

library(cli)
library(pak)
library(Giotto)
library(arrow)
library(future.apply)
library(Seurat)
#library(SeuratDisk)
library(ggplot2)
library(ggrepel)
library(dplyr)


############################################
############################################
# Preprocess combined object

gi_py <- "/your/giotto_env/bin/python"

Sys.setenv(RETICULATE_PYTHON = gi_py)  # set first
library(reticulate)
use_python(gi_py, required = TRUE)
py_config()  # sanity check

instrs <- createGiottoInstructions(python_path = gi_py) 

#load existing Giotto object
og_g <- loadGiotto("/path/to/Giotto/object/")

g <- filterGiotto(gobject = og_g,
                  spat_unit = c("aggregate", "z0", "z1", "z2", "z3", "z4", "z5", "z6"), 
                  expression_threshold = 1,
                  feat_det_in_min_cells = 5,
                  min_det_feats_per_cell = 5,
                  spat_unit_fsub = ":all:") # Cells removed 10496 out of 184911, Feats removed 0 out of 500


# normalize on aggregated layer
g <- normalizeGiotto(gobject = g, spat_unit = 'aggregate',
                     scalefactor = 1000, verbose = T)

g <- normalizeGiotto(gobject = g, spat_unit = 'aggregate',
                     norm_methods = 'pearson_resid', name = "scaled")

g <- addStatistics(gobject = g, spat_unit = 'aggregate')


#Plot aggregate total expression levels
spatPlot2D(gobject = g, spat_unit = 'aggregate',
           cell_color = 'total_expr', color_as_factor = F,
           point_size = 1, point_alpha = 0.5, coord_fix_ratio = T,
           save_param=list(save_name = file.path('Figures/SpatPlot2D_merfish_aggregate'), dpi=600),
           save_plot=F)


# Load in Single Cell Data and set up for single cell label transfer ####

# load existing seurat rds object
seurat_object <- readRDS("/projectnb/rd-spat/DATA/Collab/Murphy_lab/sc_merge_seurat.rds")

expression_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
giotto_SC_join <- createGiottoObject(expression = expression_matrix)
SO_meta <- seurat_object@meta.data
SO_meta$cell_ID <- rownames(SO_meta)
rownames(SO_meta) <- NULL

#add metadata
giotto_SC_join <- addCellMetadata(giotto_SC_join, 
                                  new_metadata = SO_meta, 
                                  by_column = TRUE, 
                                  column_cell_ID = "cell_ID")
sca <- giotto_SC_join

# match feats       
ufids <- intersect(featIDs(g), featIDs(giotto_SC_join)) # 491 feats in common

# strip gobjects down
sc_join <- giotto() |> 
  setGiotto(giotto_SC_join[[c("expression", "spatial_locs"), "raw"]]) 


#this can be improved still (ie using getExpression() function, for example, rather than explicitly pulling out from the object...)
x_join <- createGiottoObject(
  expression = g@expression[["aggregate"]][["rna"]][["raw"]]@exprMat,
  spatial_locs = g@spatial_locs[["aggregate"]][["raw"]]@coordinates
)

j <- joinGiottoObjects(
  list(sc_join[ufids], x_join[ufids]), 
  gobject_names = c("sc", "x")) 

j <- filterGiotto(j,
                  expression_threshold = 1, 
                  feat_det_in_min_cells = 1, 
                  min_det_feats_per_cell = 10)

j <- normalizeGiotto(j)  


# PCA
j <- runPCA(j) # default is normalized (not scaled); here, all genes used 
screePlot(j, ncp = 50)

# UMAP
j <- runUMAP(j, dimensions_to_use = 1:25) 


# Integration of scRNA data using Harmony
j <- runGiottoHarmony(j,
                      vars_use = "list_ID",
                      dim_reduction_name = "pca",
                      dimensions_to_use = 1:25,
                      name = "harmony_standard_25"
) #converged after 3 iterations


j <- runUMAP(j, dim_reduction_to_use = "harmony",
             dim_reduction_name = "harmony_standard_25",
             name = "hstandard_umap_25",
             dimensions_to_use = 1:25)  



# Check scRNAseq overlaid onto Merscope data
plotUMAP(j, cell_color = "list_ID", dim_reduction_name = "hstandard_umap_25", #hstandard_umap_25
         cell_color_code = c("red"), point_size = 0.3,
         point_border_stroke = 0, select_cell_groups = "sc", other_point_size = 0.3)


# Add metadata separately since there was an issue with metadata joining
sca_meta <- pDataDT(sca)
sca_meta[, cell_ID := paste0("sc-", cell_ID)]
j <- addCellMetadata(j, new_metadata = sca_meta, by_column = TRUE)


# Check scRNAseq-based annotations
plotUMAP(j, cell_color = "celltype_level2", dim_reduction_name = "hstandard_umap_25",
         select_cells = sca_meta$cell_ID, point_size = 0.4, point_border_stroke = 0,
         other_point_size = 0.3, other_cell_color ="#898989")


# Perform single cell label transfer to spatial object ####

# 'celltype_level2' can be changed to other levels of granularity from the scRNA data
j <- labelTransfer(j, 
                   source_cell_ids = spatIDs(j, subset = (list_ID == "sc")),
                   k = 10,
                   labels = "celltype_level2",
                   reduction_method = "harmony",
                   reduction_name = "harmony_standard_25", 
                   dimensions_to_use = 1:25
) 


# adding back metadata... 
ann_meta <- pDataDT(j)
ann_meta[, cell_ID := as.character(cell_ID)] 
ann_meta <- ann_meta[grepl("^x", cell_ID)]
ann_meta$cell_ID <- gsub('^x-', '', ann_meta$cell_ID)
ann_meta <- ann_meta[, .(cell_ID, trnsfr_celltype_level2, trnsfr_celltype_level2_prob)]


## get existing metadata
meta_existing <- pDataDT(g)

# merge with new annotation info
meta_merged <- merge(meta_existing, ann_meta, by = "cell_ID", all.x = TRUE)


# Add merged metadata back into the aggregate layer
g <- addCellMetadata(
  gobject = g,
  new_metadata = meta_merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  spat_unit = "aggregate"  
)

y <- g # set a checkpoint

# continue analysis with spatial object (newly annotated from labelTransfer)

meta <- pDataDT(y, spat_unit = "aggregate")

# create sanitized version of the cell type labels
meta[, trnsfr_celltype_level2_sanitized := gsub("[^A-Za-z0-9]", "_", trnsfr_celltype_level2)]


# add new column back into the object
y <- addCellMetadata( 
  gobject = y, 
  spat_unit = "aggregate",
  new_metadata = meta,
  by_column = TRUE
)

# dimension reduction
y <- runPCA(gobject = y,
            spat_unit = 'aggregate',
            expression_values = 'scaled',
            scale_unit = T, center = T, 
            feats_to_use= NULL) #use all genes

screePlot(y,
          ncp = 20,
          spat_unit = 'aggregate')

plotPCA(y,
        spat_unit = 'aggregate',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)

# UMAP
y <- runUMAP(y, spat_unit = 'aggregate', dimensions_to_use = 1:25, n_threads = 4)
plotUMAP(y, spat_unit = "aggregate", cell_color = "trnsfr_celltype_level2_sanitized", point_size = 0.3,
         point_border_stroke = 0, other_point_size = 0.3)


# Rebuild object before subsetting
cat("Rebuilding object to ensure perfect consistency before analysis...\n")

# Extract the final, clean, aggregated data
expression_final <- Giotto::getExpression(y, spat_unit = "aggregate", values = "normalized", output = "matrix")
metadata_final <- pDataDT(y, spat_unit = "aggregate")
spatlocs_final <- getSpatialLocations(y, spat_unit = "aggregate", output= "data.table")

# Build a new, clean object. This is the object to use for future steps
gobject_clean <- createGiottoObject(
  expression = expression_final,
  spatial_locs = spatlocs_final,
  cell_metadata = metadata_final
)

# Proceed with analysis using clean object
# Note: The default spat_unit is now "cell", which we use from now on


# (optional) Split the clean object to use downstream
sample1 <- pDataDT(gobject_clean) %>% filter(list_ID == "g1") %>% pull(cell_ID)

gobject_A <- subsetGiotto(
  gobject = gobject_clean,
  cell_ids = sample1
)

# Create the network on the clean subset
gobject_A <- createSpatialDelaunayNetwork(
  gobject = gobject_A,
  name = "Delaunay_network", 
  spat_unit = 'cell', # Use the new default spat_unit
  minimum_k = 2, 
  maximum_distance = "auto"
)

# Handle any NA values
gobject_A@cell_metadata$cell$rna@metaDT[is.na(trnsfr_celltype_level2_sanitized), trnsfr_celltype_level2_sanitized := "Unknown"]
gobject_A@expression[["cell"]][["rna"]][["normalized"]] <- gobject_A@expression[["cell"]][["rna"]][["raw"]]

# Run the analysis (using parallel setup)
plan(multisession)


# Code initially used for selecting Ligands / Receptors ####

#Load Panel & DB
#load table of genes from our custom panel (included in supplements)
#panel_path <- "merscope_panel_filtered_genes_removed.xlsx"
#panel_genes <- read_excel(panel_path, col_names = TRUE)[[1]]
#panel_genes <- unique(panel_genes)

#load("CellChatDB.mouse.rda")
#cellchat_interactions <- CellChatDB.mouse$interaction

#check entries for panel genes we are interested in
#cellchat_interactions[grepl('Vegf', cellchat_interactions$ligand, ignore.case = TRUE),] 



#Define our set of Ligands and Receptors
ligands <- c(
  "Vegfa", "Vegfa", "Bmp2", "Bmp4", "Ccl2", 
  "Cxcl10", "Fasl", "Col1a1", "Fgf2", "Igf1", "Ifna2"
)
receptors <- c(
  "Flt1", "Kdr", "Bmpr1a", "Bmpr2", "Ccr2", 
  "Cxcr4", "Fas", "Itga1", "Fgfr1", "Igf1r", "Ifnar1"
)

interactions <- Giotto::spatCellCellcom(
  gobject = gobject_A,
  spat_unit = "cell", # Use the new default spat_unit
  spatial_network_name = "Delaunay_network",
  cluster_column = "trnsfr_celltype_level2_sanitized",
  feat_set_1 = ligands,
  feat_set_2 = receptors,
  random_iter = 1000, #lower this for faster run speed... this step can take a while, especially if running on both samples at once
  do_parallel = TRUE
)

plan(sequential)

# Output interactions table csv to save time in the future
#write.csv(interactions, "", row.names = FALSE)

# test some figures... 
# Subset for Megakaryocyte_Platelet on the left (LIGAND)
meg_as_sender <- interactions[grepl("^Megakaryocyte_Platelet--", LR_cell_comb)]

# Subset for Megakaryocyte_Platelet on the right (RECEPTOR)
meg_as_receiver <- interactions[grepl("--Megakaryocyte_Platelet$", LR_cell_comb)]


#new cell pair list via sorting interactions table on adj-p and logfc values to try to find the most significant/informative interactions
#these cell pairs were taken from ~top50 (top56) entries in the (full, two-sample) interactions table (mostly in order from the top, but skipping a few here and there to get 
#a representative group of LR pairs)
cell_pairs <- c(
"ASMC--IM",               #Ccl2-Ccr2, p=0, logfc=0.0515
"Sox9_Epi--iMON",         #Ccl2-Ccr2, p=0, logfc=0.0428
"Mesothelial--CAP1_EPC",  #Cxcl10-Cxcr4, p=0, logfc=0.040
"Mesothelial--maDC",      #Fasl-Fas, p=0, logfc=0.035
"ILC--VSMC",              #Igf1-Igf1r, p=0, logfc=0.029
"IM--NK",                 #Igf1-Igf1r, p=0, logfc=0.027
"AT2--Ciliated",          #Fgf2-Fgfr1, p=0, logfc=0.026
"cDC2--B",                #Ifna2-Ifnar1, p=0, logfc=0.024
"AT1--Secretory",         #Fgf2-Fgfr1, p=0, logfc=0.023
"CAP1_EPC--Basal",        #Bmp4-Bmpr2, p=0, logfc=0.022
"VSMC--Neutrophil",       #Bmp4-Bmpr2, p=0, logfc=0.015
"Ciliated--AF1",          #Vegfa-Kdr, p=0, logfc=0.013
"Deuterosomal--ASMC",     #Vegfa-Flt1, p=0, logfc=0.012
"PMP--cDC2",              #Vegfa-Kdr, p=0, logfc=-0.009
"AF2--Megakaryocyte_Platelet", #Col1a1-Itga1, p=0, logfc=0.0106
"Neutrophil--Sox9_Epi",   #Bmp2-Bmpr1a, p=0.002, logfc=0.022
"Megakaryocyte_Platelet--CAP1_EPC" #Bmp4-Bmpr2, p=0.002, logfc=0.015
  )


# Subset the interactions table for matching LR_cell_comb entries
curated_interactions <- interactions[interactions$LR_cell_comb %in% cell_pairs, ]



# Plot function used for Figure 5A
plotCCcomDotplot(
  gobject_A,
  curated_interactions, #meg_as_sender, meg_as_receiver
  selected_LR = NULL,
  selected_cell_LR = NULL, #can filter here...
  show_LR_names = TRUE,
  show_cell_LR_names = TRUE,
  cluster_on = c("LR_expr"), #"PI","log2fc"
  cor_method = c("pearson"),
  aggl_method = c("ward.D"),
  dot_color_gradient = NULL,
  gradient_style = c("divergent"),
  show_plot = TRUE,
  return_plot = NULL,
  save_plot = NULL,
  save_param = list(base_width = 7, base_height = 8),
  default_save_name = "plotCCcomDotplot"
) 



# Visualize specific zoomed regions

g@instructions[["save_plot"]] <- FALSE
g@instructions[["return_plot"]] <- TRUE
g@instructions[["show_plot"]] <- TRUE


celltype_to_color <- c(
    ## Alveolar epithelial (cyans/teals)
    "AT1"         = "#1DE9B6",  # cyan
    "AT2"         = "#00E5FF",  # aqua
    "AT1/AT2"     = "#00B8D4",  # deep cyan
    "AEC"         = "#4DD0E1",  # sky blue
    
    ## Immune (reds/oranges/yellows)
    "Megakaryocyte"           = "#FFFFFF",  # white
    "Megakaryocyte/Platelet"  = "#F8F9FA",  # off-white
    "AM"          = "#FF0000",  # fire engine red
    "IM"          = "#FF4500",  # orange-red
    "iMON"        = "#FF6F00",  # dark orange
    "Neutrophil"  = "#FF8F00",  # amber-orange
    "Basophil"    = "#FFA500",  # vivid orange
    "Mast"        = "#FFB732",  # goldenrod
    "CD4 T"       = "#FFCA28",  # saffron
    "CD8 T"       = "#FFD54F",  # lemon yellow
    "Treg"        = "#FFF176",  # pastel yellow
    "NK"          = "#FFFF00",  # pure yellow
    "ILC"         = "#FEE101",  # bright golden yellow
    "B"           = "#F4C430",  # rich golden
    "Plasma"      = "#E6AC00",  # mustard
    "cDC1"        = "#E07B00",  # tangerine
    "cDC2"        = "#D95F02",  # burnt orange
    "maDC"        = "#BF360C",  # brick red
    "PMP"         = "#A52A2A",  # brown
    
    ## Capillary / endothelial (blues/purples)
    "CAP1/EPC"    = "blue",  # indigo
    "CAP2"        = "#00B0FF",  # bright blue
    "LEC"         = "#7C4DFF",  # purple
    "VEC"         = "#651FFF",  # violet
    "VSMC"        = "#9575CD",
    "Pericyte"    = "#5E35B1",
    
    ## Stromal / fibroblasts (greens)
    "SCMF"        = "#00BFA5",  # jade
    "ASMC"        = "#00C853",  # emerald
    "Mesothelial" = "#00FF00",  # neon green
    
    ## Epithelial - airway (pinks/magentas)
    "Secretory"       = "#FF4081",
    "Ciliated"        = "#F50057",
    "Deuterosomal"    = "#EC407A",
    "Sox9 Epi"        = "#FF80AB",
    "Basal"           = "#D81B60",
    "PNEC"            = "#FF69B4",
    
    ## Ambiguous
    "AF1" = "#8D6E63",
    "AF2" = "#A1887F",
    
    ## NA
    "NA" = "#BDBDBD"
)

# define region to Zoom in on 
zoom <- subsetGiottoLocs(
  g,
  spat_unit = ":all:",
  x_min = 4000,
  x_max = 6500,
  y_min = 2500,
  y_max = 5000, 
  z_min = NULL,
  z_max = NULL)

# plot specific ligands and receptors
insitu_plot <- spatInSituPlotPoints(
  zoom,
  feat_type = "rna",
  feats = list('rna'= c("Bmp4", "Bmpr2")),
  feats_color_code = c("Bmp4"="yellow", "Bmpr2"="green"),
  point_size = 0.35,
  polygon_feat_type = "aggregate",
  polygon_fill = "trnsfr_celltype_level2", 
  polygon_fill_as_factor = TRUE,
  polygon_fill_code = celltype_to_color,
  polygon_line_size = 0, 
  polygon_alpha = 0.5,
  save_param = list(base_height = 10, #adjust the height and width for your sample
                    base_width = 16)) 