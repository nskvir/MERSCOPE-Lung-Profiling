#Giotto analysis script


#pak::pkg_install("drieslab/Giotto")

library(cli)
library(pak)
library(Giotto)
#library(bit64)
library(arrow)
library(future.apply)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggrepel)



#Specify root directory & file paths for the two sample folders
root_dir <- "/directory/containing/sample/data/"

runs <- list(
        A = file.path(root_dir, "202403242054_VA00340-BUMC-CellBoundary-GroupA_VMSC01801", "region_0"),
        B = file.path(root_dir, "202403242055_VA00340-BUMC-CellBoundary-GroupB_VMSC03501", "region_0")
    )

# Set seed to ensure reproducibility of figures with stochastic elements
set.seed(123)


##################################################
### 1.) Load in the data, create Giotto objects
##################################################
setwd(dirname(runs$A))

##### First, create Giotto object from sample A data

#Define paths to subdirectories within the data folder
im_subdir <- file.path(runs$A, "images")
tx_path <- file.path(runs$A, "detected_transcripts.csv")
poly_path <- file.path(runs$A, "cell_boundaries.parquet")

# transcripts
tx <- data.table::fread(tx_path)
gpoints = createGiottoPoints(
    tx,
    x_colname = 'global_x',
    y_colname = 'global_y',
    feat_ID_colname = 'gene',
    split_keyword = list("Blank"),
    feat_type = c("rna", "blank")
)

# polys
polygons = readPolygonVizgenParquet(
    file = poly_path,
    calc_centroids = TRUE
)

# images
im_name_fmt <- "mosaic_%s_z%d.tif"
im_type <- list(
    pt = "PolyT",
    cb1 = "Cellbound1",
    cb2 = "Cellbound2",
    cb3 = "Cellbound3",
    dp = "DAPI"
)

# load images
imgs <- list()
for (type in (im_type)) {
    z <- 0:6
    imgs[[type]] <- createMerscopeLargeImage(
        file.path(im_subdir, sprintf(im_name_fmt, type, z)),
        transforms_file = file.path(im_subdir, "micron_to_mosaic_pixel_transform.csv"),
        name = sprintf("%s_z%d", type, z)
    )
    names(imgs[[type]]) <- sprintf("z%d", z)
}

# flip vector data ############
img_ext <- ext(imgs$PolyT$z0)

# find image y midline from first image
y_mid <- mean(c(img_ext[3], img_ext[4]))

gpoints <- lapply(
    gpoints,
    flip,
    y0 = y_mid
)
polygons <- lapply(
    polygons,
    flip,
    y0 = y_mid
)

#create giotto object
g1 <- giotto()
g1 <- setGiotto(g1, polygons)
g1 <- setGiotto(g1, gpoints)
g1 <- setGiotto(g1, imgs)

layers <- names(g1@spatial_info)
layers <- layers[layers != 'aggregate']
layers <- layers[grepl("z[0-9]", layers)] # Keeps only the Z-stacks

for (layer in seq(0, length(layers)-1)){
  g1 = calculateOverlap(g1,
                                spatial_info = paste0('z', layer),
                                feat_info = 'rna',
                                feat_subset_column = 'global_z',
                                feat_subset_ids = layer)
  
  g1 = overlapToMatrix(g1,
                               poly_info = paste0('z', layer),
                               feat_info = 'rna',
                               name = 'raw')
}

# aggregate information from multiple polygons in different z-stacks
# into a single average location/cell
g1 = aggregateStacks(gobject = g1,
                           spat_units = layers,
                           feat_type = 'rna',
                           values = 'raw',
                           summarize_expression = 'sum',
                           summarize_locations = 'mean',
                           new_spat_unit = 'aggregate')



##### Next, create a second Giotto object from sample B data
setwd(dirname(runs$B))

im_subdir <- file.path(runs$B, "images")
tx_path <- file.path(runs$B, "detected_transcripts.csv")
poly_path <- file.path(runs$B, "cell_boundaries.parquet")

# transcripts
tx <- data.table::fread(tx_path)
gpoints = createGiottoPoints(
    tx,
    x_colname = 'global_x',
    y_colname = 'global_y',
    feat_ID_colname = 'gene',
    split_keyword = list("blank"),
    feat_type = c("rna", "blank")
)

# polys
polygons = readPolygonVizgenParquet(
    file = poly_path,
    calc_centroids = TRUE
)

# images
im_name_fmt <- "mosaic_%s_z%d.tif"
im_type <- list(
    pt = "PolyT",
    cb1 = "Cellbound1",
    cb2 = "Cellbound2",
    cb3 = "Cellbound3",
    dp = "DAPI"
)

# load images
imgs <- list()
for (type in (im_type)) {
    z <- 0:6
    imgs[[type]] <- createMerscopeLargeImage(
        file.path(im_subdir, sprintf(im_name_fmt, type, z)),
        transforms_file = file.path(im_subdir, "micron_to_mosaic_pixel_transform.csv"),
        name = sprintf("%s_z%d", type, z)
    )
    names(imgs[[type]]) <- sprintf("z%d", z)
}

# flip vector data
img_ext <- ext(imgs$PolyT$z0)

# find image y midline from first image
y_mid <- mean(c(img_ext[3], img_ext[4]))

gpoints <- lapply(
    gpoints,
    flip,
    y0 = y_mid
)
polygons <- lapply(
    polygons,
    flip,
    y0 = y_mid
)



# create Giotto object
g2 <- giotto()
g2 <- setGiotto(g2, polygons)
g2 <- setGiotto(g2, gpoints)
g2 <- setGiotto(g2, imgs)

layers <- names(g2@spatial_info)
layers <- layers[layers != 'aggregate']
layers <- layers[grepl("z[0-9]", layers)] # Keeps only the Z-stacks

for (layer in seq(0, length(layers)-1)){
  g2 = calculateOverlap(g2,
                                spatial_info = paste0('z', layer),
                                feat_info = 'rna',
                                feat_subset_column = 'global_z',
                                feat_subset_ids = layer)
  
  g2 = overlapToMatrix(g2,
                               poly_info = paste0('z', layer),
                               feat_info = 'rna',
                               name = 'raw')
}

# aggregate information from multiple polygons in different z-stacks
# into a single average location/cell
g2 = aggregateStacks(gobject = g2,
                           spat_units = layers,
                           feat_type = 'rna',
                           values = 'raw',
                           summarize_expression = 'sum',
                           summarize_locations = 'mean',
                           new_spat_unit = 'aggregate')



# Combine the two spatial objects (can be memory intensive)
g <- joinGiottoObjects(c(g1,g2), gobject_names=c("g1","g2"), x_shift=c(0,10000), y_shift=c(0,0))



######################################
### 2.) Filtering and normalization
######################################

g <- filterGiotto(gobject = g,
                         spat_unit = 'aggregate',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5,
                         poly_info = c('aggregate'))

# Number of cells removed:  10449  out of  184911 
# Number of feats removed:  0  out of  550



# normalize on aggregated layer
g <- normalizeGiotto(gobject = g, spat_unit = 'aggregate',
                            scalefactor = 1000, verbose = T)
g <- addStatistics(gobject = g, spat_unit = 'aggregate')
g <- normalizeGiotto(gobject = g, spat_unit = 'aggregate',
                            norm_methods = 'pearson_resid', name = "scaled")


#Plot aggregate total expression levels with background image
spatPlot2D(gobject = g, spat_unit = 'aggregate',
           cell_color = 'total_expr', color_as_factor = F,
           image_name = 'g1-DAPI_z1', show_image = FALSE,
           point_size = 1, point_alpha = 0.5, coord_fix_ratio = T,
           save_param=list(save_name = file.path('Figures/SpatPlot2D_merfish_aggregate'), dpi=600),
           save_plot=F)



################################################################################################
### 3.) Merge polygons from identified megakaryocytes (identified with hex-binning)
################################################################################################


### NOTE: The block of commented code below shows how we obtained the cell IDs of polygons comprising 
### our MKs, followed by how we initially merged them and redrew boundaries around them. This process is streamlined
### in the uncommented code (all IDs are already collected, and boundaries are loaded in via .shp files for reproducibility.)

### Each of the MK coordinates and cell IDs below was obtained through manual inspection of candidate regions,
### based on hex-binning in our previous script. 


#####
### Example code illustrating how to view cell IDs in a given region (MK1 in this case):
#####

#cell_loc_data = combineCellData(gobject = g, include_poly_info = F, poly_info = 'aggregate')
#cell_loc_data = cell_loc_data$rna
#subset_cell_loc_data = cell_loc_data[sdimx > 4400 & sdimx < 4700 & sdimy > 6200 & sdimy < 6400] #GroupA, FOV2
#poly_data = combineCellData(gobject = g, include_poly_info = T, poly_info = 'aggregate')
#poly_data = poly_data$rna
#subset_poly_data = poly_data[sdimx > 4400 & sdimx < 4700 & sdimy > 6200 & sdimy < 6400] #GroupA, FOV2
#pl = ggplot()
#pl = pl + geom_polygon(data = subset_poly_data, aes(x = x, y = y, fill = 'red', group = cell_ID))
#pl = pl + geom_point(data = subset_cell_loc_data, aes(x = sdimx, y = sdimy, size = 3))
#pl = pl + geom_text_repel(data = subset_cell_loc_data, aes(x = sdimx, y = sdimy, label = cell_ID, size = 3))
#pl
#####

# Sample A coordinates:

# MK1 (Group A, FOV2), coordinates: x(4400,4700), y(6200,6400)
#cell_ids_to_merge = c('g1-1864822800107100045', 
                      #'g1-1864822800107100060',
                      #'g1-1864822800107100067', 
                      #'g1-1864822800107200047', 
                      #'g1-1864822800107100070', 
                      #'g1-1864822800107100088')

# MK2 (Group A, FOV3), coordinates: x(1800,2100), y(5000,5200)
#cell_ids_to_merge = c('g1-1864822800124200214', 
                      #'g1-1864822800124200237', 
                      #'g1-1864822800124100282', 
                      #'g1-1864822800124100271')

# MK3 (Group A, FOV4), coordinates: x(2950,3250), y(1750,1950)
#cell_ids_to_merge = c('g1-1864822800193200301', 
                      #'g1-1864822800193200304', 
                      #'g1-1864822800193100449', 
                      #'g1-1864822800193100457')

# MK4 (Group A, FOV5), coordinates: x(3900,4250), y(1300,1500)
#cell_ids_to_merge = c('g1-1864822800205100469', 
                      #'g1-1864822800205100464')


# Sample B coordinates:

# MK5 (Group B, FOV2), coordinates (x-shifted +10k from merged samples): x(15850,16050), y(11700,11825)
#cell_ids_to_merge = c('g2-1865258400051100139', 
                      #'g2-1865258400051100140')

# MK6 (Group B, FOV4), coordinates (x-shifted +10k from merged samples): x(14200,14400), y(10050,10250)
#cell_ids_to_merge = c('g2-1865258400108100563', 
                      #'g2-1865258400108100559')

# MK7 (Group B, FOV6), coordinates (x-shifted +10k from merged samples): x(15980,16100), y(9500,9650)
#cell_ids_to_merge = c('g2-1865258400131100564', 
                      #'g2-1865258400131100575',
                      #'g2-1865258400131200382')

# MK8 (Group B, FOV31), coordinates (x-shifted +10k from merged samples): x(17350,17500), y(2700,2800)
#cell_ids_to_merge = c('g2-1865258400394100872', 
                      #'g2-1865258400394200535')


### For each set of cell IDs above, we can draw a new boundary around them, using Terra (as below)

#test = getPolygonInfo(g, polygon_name = 'aggregate')
#original = test[test$poly_ID %in% cell_ids_to_merge]
#plot(original)
#MK1 = terra::draw(x = 'polygon')
#MK1_centroid = terra::centroids(MK1)
#plot(MK1)
#plot(MK1_centroid, add = TRUE)

#MK1 = terra::makeValid(MK1)

#terra::values(MK1) = data.table::data.table(poly_ID = 'MK_1',
                                            #mean_z_level = 6,
                                            #stack = NA, agg_n = 7, 
                                            #valid = TRUE)

### Lastly, the new boundary can be exported as a .shp file and loaded in later

#writeVector(MK1, "groupA_mk1.shp", overwrite=FALSE)             


#ALL Meg IDS combined to merge:
cell_ids_to_merge = c('g1-1864822800107100045', 
                      'g1-1864822800107100060', 
                      'g1-1864822800107100067', 
                      'g1-1864822800107200047', 
                      'g1-1864822800107100070', 
                      'g1-1864822800107100088',
                      'g1-1864822800124200214', 
                      'g1-1864822800124200237', 
                      'g1-1864822800124100282', 
                      'g1-1864822800124100271',
                      'g1-1864822800193200301', 
                      'g1-1864822800193200304', 
                      'g1-1864822800193100449', 
                      'g1-1864822800193100457',
                      'g1-1864822800205100469', 
                      'g1-1864822800205100464',
                      'g2-1865258400051100139',
                      'g2-1865258400051100140',
                      'g2-1865258400108100563',
                      'g2-1865258400108100559',
                      'g2-1865258400131100564',
                      'g2-1865258400131100575',
                      'g2-1865258400131200382',
                      'g2-1865258400394100872',
                      'g2-1865258400394200535') #all of them together, to use to generate 'remaining'

#show original plot with locations of cells to be merged
test = getPolygonInfo(g, polygon_name = 'aggregate')
original = test[test$poly_ID %in% cell_ids_to_merge]
plot(original)


library(terra)

#to read drawn vectors back in:
MK1 <- vect("groupA_mk1.shp")
MK2 <- vect("groupA_mk2.shp")
MK3 <- vect("groupA_mk3.shp")
MK4 <- vect("groupA_mk4.shp")
MK5 <- vect("groupB_mk5.shp")
MK6 <- vect("groupB_mk6.shp")
MK7 <- vect("groupB_mk7.shp")
MK8 <- vect("groupB_mk8.shp")

MK1_centroid <- terra::centroids(MK1)
MK2_centroid <- terra::centroids(MK2)
MK3_centroid <- terra::centroids(MK3)
MK4_centroid <- terra::centroids(MK4)
MK5_centroid <- terra::centroids(MK5)
MK6_centroid <- terra::centroids(MK6)
MK7_centroid <- terra::centroids(MK7)
MK8_centroid <- terra::centroids(MK8)


# get other polygons
remaining = test[!test$poly_ID %in% cell_ids_to_merge] 

# combine other polygons with new MK polygons
new_polygons = rbind(remaining, MK1, MK2, MK3, MK4, MK5, MK6, MK7, MK8)

# create new giotto polygon object that can be added to the giotto object
new_polygons_giotto = createGiottoPolygon(x = new_polygons, name = 'newaggregate', calc_centroids = TRUE)

g <- setGiotto(g, new_polygons_giotto)

# create spatial locations
cell_info = terra::values(g@spatial_info$newaggregate@spatVectorCentroids)
xydims = terra::geom(g@spatial_info$newaggregate@spatVectorCentroids)
new_coords = data.table::data.table(cell_ID = cell_info$poly_ID, sdimx = xydims[,3], sdimy = xydims[,4])
new_locations_giotto = createSpatLocsObj(coordinates = new_coords, name = 'raw', spat_unit = 'newaggregate')

g <- setGiotto(g, new_locations_giotto)

showGiottoSpatialInfo(g)

pDataDT(g, 'newaggregate')

plot = spatInSituPlotPoints(g,
                            show_polygon = TRUE,
                            spat_unit = 'newaggregate',
                            polygon_feat_type = 'newaggregate',
                            polygon_color = 'white',
                            polygon_line_size = 0.1,
                            polygon_alpha = 1,
                            polygon_fill = NULL,
                            polygon_fill_as_factor = T,
                            coord_fix_ratio = T, 
                            return_plot = T)

#check if MK1 is merged
#plot = plot + ggplot2::xlim(c(4400, 4700)) + ggplot2::ylim(c(6200, 6400))
#plot


#recalculate transcript overlaps with newaggregate
g = calculateOverlap(g,
                                spatial_info = 'newaggregate',
                                feat_info = 'rna')
  
g = overlapToMatrix(g,
                               poly_info = 'newaggregate',
                               feat_info = 'rna',
                               name = 'raw')

#filter again
g <- filterGiotto(gobject = g,
                         spat_unit = 'newaggregate',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5,
                         poly_info = c('newaggregate'))

# normalize on newaggregate layer
g <- normalizeGiotto(gobject = g, spat_unit = 'newaggregate',
                            scalefactor = 1000, verbose = T)
g <- addStatistics(gobject = g, spat_unit = 'newaggregate')
g <- normalizeGiotto(gobject = g, spat_unit = 'newaggregate',
                            norm_methods = 'pearson_resid', update_slot = 'pearson')


##########################################################
### 4.) Load single cell reference data from LungMAP
##########################################################


# The commented steps were used to convert the provided data into a usable format

# Single-cell data was obtained from the LungMAP repository here:
# https://www.lungmap.net/dataset/?dataset_id=LMEX0000004397

# To convert the provided .h5ad files to .h5seurat format, the following repository, 'seurat-disk' was used:
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

# The following commented steps are for reference
# They can be skipped if simply loading in the .h5seurat file provided

### Converting .h5ad file into an .h5Seurat file...
#Convert("/path/to/LungMAP_MouseLung_CellRef.v1.1.h5ad", dest = "./data.h5seurat", overwrite = TRUE)

### Read in the corresponding metadata (this gives celltype labels, etc)
#cell_metadata <- read.table("/path/to/Cell.Metadata.txt", header=TRUE)
#seurat_object <- AddMetaData(seurat_object, metadata=cell_metadata)

### Read the UMAP coordinates file 
#umap_coordinates <- read.table("/path/to/UMAP-coordinates.txt", header=FALSE, row.names= 1)

### Add UMAP embeddings to the Seurat object
###NOTE:  Ensure cell names/barcodes match: The cell names (row names in the UMAP file and metadata) must match the cell names in the Seurat object
#seurat_object[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap_coordinates), key = "UMAP_")


#load the converted (.h5Seurat) file into Seurat
seurat_object <- LoadH5Seurat("/location/of/data.h5seurat")

expression_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
giotto_SC_join <- createGiottoObject(expression = expression_matrix)
SO_meta <- seurat_object@meta.data

#check order of cell_IDs
giotto_cell_IDs <- giotto_SC_join@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]]
seurat_cell_IDs <- rownames(SO_meta)
all_match <- all(giotto_cell_IDs == seurat_cell_IDs) #True...
mismatches <- which(giotto_cell_IDs != seurat_cell_IDs) #empty...

giotto_SC_join <- addCellMetadata(giotto_SC_join, 
                                  new_metadata = SO_meta)
sca <- giotto_SC_join

# Match feats between 500-gene panel and scRNA reference   
ufids <- intersect(featIDs(g), featIDs(giotto_SC_join)) # 491 feats in common

# strip gobjects down
sc_join <- giotto() |> 
  setGiotto(giotto_SC_join[[c("expression", "spatial_locs"), "raw"]]) 


x_join <- createGiottoObject(
  expression = g@expression[["newaggregate"]][["rna"]][["raw"]]@exprMat,
  spatial_locs = g@spatial_locs[["newaggregate"]][["raw"]]@coordinates
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
j <- runPCA(j) #all genes used 
# scree
screePlot(j, ncp = 50) 

# UMAP
j <- runUMAP(j, dimensions_to_use = 1:25) 

# harmony runs
j <- runGiottoHarmony(j,
                      vars_use = "list_ID",
                      dim_reduction_name = "pca",
                      dimensions_to_use = 1:25,
                      name = "harmony_standard_25"
    )

j <- runUMAP(j, dim_reduction_to_use = "harmony",
             dim_reduction_name = "harmony_standard_25",
             name = "hstandard_umap_25",
             dimensions_to_use = 1:25)  



# scRNAseq vs Merscope data
plotUMAP(j, cell_color = "list_ID", dim_reduction_name = "hstandard_umap_25", #hstandard_umap_25
         cell_color_code = c("red"), point_size = 0.3,
         point_border_stroke = 0, select_cell_groups = "sc", other_point_size = 0.3)

sca_meta <- pDataDT(sca)
sca_meta[, cell_ID := paste0("sc-", cell_ID)]
j <- addCellMetadata(j, new_metadata = sca_meta, by_column = TRUE)


# Draw scRNAseq annotations vs MERSCOPE data (Figure 2A)
plotUMAP(j, cell_color = "celltype_level2", dim_reduction_name = "hstandard_umap_25",
         select_cells = sca_meta$cell_ID, point_size = 0.4, point_border_stroke = 0,
         other_point_size = 0.3, other_cell_color ="#898989")


#check distribution of cells and feats per cell
#filterDistributions(j)
#filterDistributions(j, detection = 'cells')


#########################
### 5.) Label transfer
#########################


#'celltype_level2' can be changed to whichever level of granularity we want from the reference data
j <- labelTransfer(j, 
                   source_cell_ids = spatIDs(j, subset = (list_ID == "sc")),
                   k = 10,
                   labels = "celltype_level2",
                   reduction_method = "harmony",
                   reduction_name = "harmony_standard_25", 
                   dimensions_to_use = 1:25
) 

#Check transfer probability distribution
#hist(pDataDT(j)$trnsfr_celltype_level2_prob,
     #breaks = 50,
     #col = "skyblue",
     #main = "Transfer Probability Distribution",
     #xlab = "Transfer Probability",
     #ylab = "Number of Cells")

# Combined clusters after label transfer (Figure 2B)
plotUMAP(j, 
         cell_color = "trnsfr_celltype_level2", #trnsfr_cell_types (according to how labeled)
         dim_reduction_name = "hstandard_umap_25", 
         #cell_color_code = color_mapping, #change value for number of clusters 
         point_size = 1, 
         # select_cells = sca_meta$cell_ID,
         # other_point_size = 0.3,
         point_border_stroke = 0
) 

# adding back metadata
ann_meta <- pDataDT(j)
ann_meta[, cell_ID := as.character(cell_ID)] 
ann_meta <- ann_meta[grepl("^x", cell_ID)]
ann_meta$cell_ID <- gsub('^x-', '', ann_meta$cell_ID)
ann_meta <- ann_meta[, .(cell_ID, trnsfr_celltype_level2, trnsfr_celltype_level2_prob)]

# get existing metadata from newaggregate
meta_existing <- pDataDT(g, spat_unit = "newaggregate")

# merge with new annotation info
meta_merged <- merge(meta_existing, ann_meta, by = "cell_ID", all.x = TRUE)

# Add merged metadata back into the 'newaggregate' layer
g <- addCellMetadata(
  gobject = g,
  new_metadata = meta_merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  spat_unit = "newaggregate"  # specify newaggregate
)

#save with new name
y <- g



#colors <- getDistinctColors(39)


#color scheme used in Fig.1
celltype_to_color <- c(
  ## Alveolar epithelial (cyan)
  "AT1"              = "#4ee7b6", # Alveolar type I
  "AT2"              = "#4fe2fd", # Alveolar type II
  "AT1/AT2"          = "#3fb6d2", # Transitional alveolar state
  
  ## Capillary / endothelial (blues, violet)
  "AEC"              = "#65cedf", # Arterial endothelial cell
  "CAP1/EPC"         = "#3000fc", # 
  "CAP2"             = "#43adfd", #
  "LEC"              = "#804cfc", # Lymphatic endothelial cell
  "VEC"              = "#6b1afc", # Venous endothelial cell
    
  ## Immune (reds, oranges, yellows)
  "Megakaryocyte/Platelet" = "#f8f8f9",  # MKs are white for highlighting
  "AM"               = "#f4261c", # Alveolar macrophage
  "IM"               = "#f54f1d", # Interstitial macrophage
  "IMON"             = "#de2d26", # Inflammatory macrophage
  "Neutrophil"       = "#f69322", #
  "Basophil"         = "#f7a723", #
  "Mast"             = "#f8b93e", #
  "CD4 T"            = "#f9cb39", #
  "CD8 T"            = "#f9d558", #
  "Treg"             = "#fbf17c", #
  "NK"               = "#fcff2d", #
  "ILC"              = "#f8e12a", # Innate lymphoid cell
  "B"                = "#eec63d", #
  "cDC1"             = "#d77d1c", # Classical dendritic cell subset 1 (myeloid APC; cross-present to CD8 T cells, secrete IL-12, support Th1/cytotoxic)
  "cDC2"             = "#d1641a", # Classical dendritic cell subset 2 (myeloid APC; present to CD4 T cells, support Th1/Th2/Th17; superior antigen presentation)
  "maDC"             = "#b63d19", # Mature dendritic cell
  

  ## Mesenchymal (greens, indigo, browns)
  "SCMF"             = "#3bbca4", # Secondary crest myofibroblast
  "ASMC"             = "#37c655", # Airway smooth muscle
  "Mesothelial"      = "#47fd21", # 
  "VSMC"             = "#9475cb", # Vascular smooth muscle
  "Pericyte"         = "#6034ae", # 
  "AF1"              = "#8a6e63", # Alveolar fibroblast 1
  "AF2"              = "#9e887f", # Alveolar fibroblast 2
  "PMP"              = "#9d302e", # Proliferative mesenchymal progenitor

  ## Epithelial airway (pinks)
  "Secretory"        = "#f54a83", #
  "Ciliated"         = "#eb2459", #
  "Deuterosomal"     = "#e2497a", #
  "Sox9 Epi"         = "#f784ac", # Sox9+ airway progenitor
  "Basal"            = "#cf2b61", #
  "PNEC"             = "#f76eb4", # Pulmonary neuroendocrine cell

  ## NA or unlabeled (gray)
  "<NA>"               = "#7e7e7e" #
)



# All clusters (Figure 1)
spatInSituPlotPoints(y,
                     show_polygon = TRUE,
                     spat_unit = 'newaggregate',
                     polygon_feat_type = 'newaggregate',
                     polygon_color = 'black',
                     polygon_line_size = 0.01,
                     polygon_fill_code = celltype_to_color,
                     polygon_alpha = 1,
                     polygon_fill = 'trnsfr_celltype_level2',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = T,
                     save_plot=F,
                     save_param = list(save_name = file.path('Figures/Leiden_Mapping'),
                                       dpi = 600))


#############################################################
### 6.) Preliminary Spatial Network & Proximity Enrichment
#############################################################


# highly variable genes
y <- calculateHVF(gobject = y,
                         spat_unit = 'newaggregate',
                         method = 'var_p_resid',
                         expression_values = 'scaled',
                         show_plot = T)

# dimension reduction
y <- runPCA(gobject = y,
                   spat_unit = 'newaggregate',
                   expression_values = 'scaled',
                   scale_unit = T, center = T, 
                   feats_to_use=NULL) #use all genes

plotPCA(y,
        spat_unit = 'newaggregate',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)


y <- runUMAP(y, dimensions_to_use = 1:25, n_threads = 4, spat_unit = 'newaggregate')

y <- addStatistics(gobject = y, spat_unit = 'newaggregate')


### create spatial network
y <- createSpatialNetwork(gobject = y, 
                                     spat_unit = 'newaggregate',
                                     minimum_k = 2, 
                                     maximum_distance_delaunay = 30)


cell_proximities = cellProximityEnrichment(gobject = y,
                                           spat_unit = 'newaggregate',
                                           cluster_column = 'trnsfr_celltype_level2',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)


### Raw output of Fig.4A
cellProximityBarplot(gobject = y, 
                     CPscore = cell_proximities, 
                     min_orig_ints = 5, 
                     min_sim_ints = 5, 
                     p_val = 0.5,
                     save_plot = F,
                     save_param = list(base_width = 9, base_height = 9,
                                       save_name = file.path('Figures/proximity_barplot_FULL'),
                                       dpi = 600)) 

# Supplemental figure heatmap
cellProximityHeatmap(gobject = y, 
                     CPscore = cell_proximities, 
                     order_cell_types = T, 
                     scale = T,
                     color_breaks = c(-1.5, 0, 1.5), 
                     color_names = c('blue', 'white', 'red'),
                     save_plot = F,
                     save_param = list(base_width = 7, base_height = 7,
                                       save_name = file.path('Figures/proximity_heatmap_FULL'),
                                       dpi = 600)) 

# Supplemental figure network
cellProximityNetwork(gobject = y, 
                     CPscore = cell_proximities, 
                     remove_self_edges = F, 
                     only_show_enrichment_edges = T,
                     save_plot = F,
                     save_param = list(base_width = 7, base_height = 7,
                                       save_name = file.path('Figures/proximity_network_FULL'),
                                       dpi = 600))




#########################################
### 7.) MK Niche - Stratified Barplots
#########################################

library(Giotto)
library(data.table)
library(ggplot2)


# Relabel our 8 curated MKs in the metadata
mk_ids <- c("MK_1","MK_2","MK_3","MK_4","MK_5","MK_6","MK_7","MK_8")

meta <- pDataDT(y, spat_unit = "newaggregate")
meta[cell_ID %in% mk_ids, trnsfr_celltype_level2 := "Megakaryocyte"]

y_renamed <- addCellMetadata(
  gobject = y,
  new_metadata = meta,
  spat_unit = "newaggregate",
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)

### Ensure a KNN spatial network exists
net_name <- "knn_mk_network"
has_net <- ("newaggregate" %in% names(y_renamed@spatial_network)) &&
           (net_name %in% names(y_renamed@spatial_network$newaggregate))
if (!has_net) {
  y_renamed <- createSpatialKNNnetwork(
    gobject = y_renamed,
    name = net_name,
    spat_unit = "newaggregate",
    k = 100,
    maximum_distance = 50,
    minimum_k = 0,
    return_gobject = TRUE
  )
}

### Pull the network and attach from/to level2 types
knn_mk_network <- getSpatialNetwork(
  y_renamed, name = net_name, output = "networkDT", spat_unit = "newaggregate"
)
meta2 <- pDataDT(y_renamed, spat_unit = "newaggregate")[, .(cell_ID, trnsfr_celltype_level2)]

# add 'from' types
meta_from <- copy(meta2); setnames(meta_from, "cell_ID", "from")
knn_mk_network[meta_from, on = .(from), trnsfr_celltype_level2_from := i.trnsfr_celltype_level2]

# add 'to' types
meta_to <- copy(meta2); setnames(meta_to, "cell_ID", "to")
knn_mk_network[meta_to, on = .(to), trnsfr_celltype_level2_to := i.trnsfr_celltype_level2]

### keep only edges that originate from MKs; bin distances
mk_neighbors <- knn_mk_network[trnsfr_celltype_level2_from == "Megakaryocyte"]

mk_neighbors[, distance_bin := cut(
  distance,
  breaks = seq(0, 50, by = 10),
  include.lowest = TRUE,
  labels = c("0–10","10–20","20–30","30–40","40–50")
)]
mk_neighbors <- mk_neighbors[!is.na(distance_bin)]

### counts -> proportions per bin (by neighbor type)
type2_counts <- mk_neighbors[, .N, by = .(distance_bin, trnsfr_celltype_level2_to)]
type2_counts[, total := sum(N), by = distance_bin]
type2_counts[, proportion := N / total]




level2_colors <- c(
  # Alveolar epithelium (greens)
  "AT1" = "#369d77", #
  "AT2" = "#419f31", #
  "AT1/AT2" = "#71c1a5", #

  # Fibroblasts & Mesenchymal cells (yellows/oranges)
  "AF1" = "#fada3f", #
  "AF2" = "#fae18f", #
  "SCMF" = "#f6b066",#
  "PMP" = "#eb7249", #

  # Endothelial (blues)
  "CAP1/EPC" = "#3777b2", #
  "CAP2" = "#abcde2",     #
  "AEC" = "#76add5",      #
  "VSMC" = "#3770b3",     #
  "Pericyte" = "#172f6a", #

  # Immune - Macrophage/Monocyte/DC (reds)
  "AM" = "#da2c27",   #
  "IM" = "#f27050",   #
  "iMON" = "#d5392e", #
  "maDC" = "#f49575", #
  "cDC1" = "#f6bda3", #

  # Immune - Lymphoid (purples/pinks)
  "CD4 T" = "#9550a2", #
  "NK" = "#b981bc",    #
  "ILC" = "#d1641a",   #
  "B" = "#de368b",     #

  # Other immune (neutrals)
  "Neutrophil" = "#9e1d1c", #
  "Basophil" = "#f68420",   #
  "Mast" = "#f8b955",       #

  # Epithelial (teal/magenta)
  "Sox9 Epi" = "#95d2c7",   #
  "Mesothelial" = "#b6e630",#

  # Platelets/MKs
  "Megakaryocyte/Platelet" = "#fefdd1", #

  # Uncategorized
  "NA" = "#7f7f7f" #
)

### harmonize labels to match palette keys; lock factor order
type2_counts[, trnsfr_celltype_level2_to := trimws(as.character(trnsfr_celltype_level2_to))]
type2_counts[is.na(trnsfr_celltype_level2_to) | trnsfr_celltype_level2_to == "<NA>",
             trnsfr_celltype_level2_to := "NA"]
type2_counts[trnsfr_celltype_level2_to == "IMON", trnsfr_celltype_level2_to := "iMON"]

type2_counts[, trnsfr_celltype_level2_to :=
               factor(trnsfr_celltype_level2_to, levels = names(level2_colors))]

## Plot (Figure 4C)
ggplot(type2_counts, aes(x = distance_bin, y = proportion, fill = trnsfr_celltype_level2_to)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  scale_fill_manual(values = level2_colors,
                    limits = names(level2_colors),
                    drop   = FALSE,
                    na.value = "#7f7f7f") +
  labs(
    title = "Specific Cell Type Proportions near Megakaryocytes",
    x = "Distance (µm)", y = "Cell Proportion (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### Now with broad celltypes (Figure 4B)

### Define fine -> broad mapping

level2_to_broad <- list(
  "AEC" = "Endothelial",
  "CAP1/EPC" = "Endothelial",
  "CAP2" = "Endothelial",
  "LEC" = "Endothelial",
  "VEC" = "Endothelial",

  "Pericyte" = "Mesenchymal",
  "VSMC" = "Mesenchymal",
  "SCMF" = "Mesenchymal",
  "ASMC" = "Mesenchymal",
  "Mesothelial" = "Mesenchymal",
  "PMP" = "Mesenchymal",

  "AF1" = "Fibroblast",
  "AF2" = "Fibroblast",

  "AT1" = "Alveolar Epithelium",
  "AT2" = "Alveolar Epithelium",
  "AT1/AT2" = "Alveolar Epithelium",

  "Sox9 Epi" = "Airway Epithelium",
  "Secretory" = "Airway Epithelium",
  "Ciliated" = "Airway Epithelium",
  "Deuterosomal" = "Airway Epithelium",
  "Basal" = "Airway Epithelium",
  "PNEC" = "Airway Epithelium",

  "IM" = "Myeloid",
  "AM" = "Myeloid",
  "iMON" = "Myeloid",
  "cDC1" = "Myeloid",
  "maDC" = "Myeloid",

  "Mast" = "Granulocyte",
  "Neutrophil" = "Granulocyte",
  "Basophil" = "Granulocyte",

  "B" = "Lymphocyte",
  "CD4 T" = "Lymphocyte",
  "NK" = "Lymphocyte",
  "ILC" = "Lymphocyte",

  "Megakaryocyte/Platelet" = "Megakaryocyte",

  "NA" = "Unlabeled/Other"
)


if (is.list(level2_to_broad)) {
  level2_to_broad <- unlist(level2_to_broad, use.names = TRUE)
}

## normalize neighbor labels to match keys, then map to broad
mk_neighbors[, trnsfr_celltype_level2_to := trimws(as.character(trnsfr_celltype_level2_to))]
mk_neighbors[is.na(trnsfr_celltype_level2_to) | trnsfr_celltype_level2_to == "<NA>", trnsfr_celltype_level2_to := "NA"]
mk_neighbors[trnsfr_celltype_level2_to == "IMON", trnsfr_celltype_level2_to := "iMON"]

mk_neighbors[, broad_category := unname(level2_to_broad[trnsfr_celltype_level2_to])]
mk_neighbors[is.na(broad_category), broad_category := "Unlabeled/Other"]

## 3) Proportions within each distance bin
broad_counts <- mk_neighbors[, .N, by = .(distance_bin, broad_category)]
broad_counts[, total := sum(N), by = distance_bin]
broad_counts[, proportion := N / total]

### broad palette + order
broad_palette <- c(
  "Alveolar Epithelium"="#00B8D4",
  "Airway Epithelium"  ="#F06292",
  "Endothelial"        ="#3F51B5",
  "Mesenchymal"        ="#00C853",
  "Fibroblast"         ="#B2FF59",
  "Myeloid"            ="#FF7043",
  "Granulocyte"        ="#FFA726",
  "Lymphocyte"         ="#9575CD",
  "Megakaryocyte"      ="#FFFDD0",
  "Unlabeled/Other"    ="#D3D3D3"
)
broad_counts[, broad_category := factor(broad_category, levels = names(broad_palette))]

### plot
ggplot(broad_counts, aes(x = distance_bin, y = proportion, fill = broad_category)) +
  geom_col() +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%")) +
  scale_fill_manual(values = broad_palette, limits = names(broad_palette), drop = FALSE) +
  labs(
    title = "Broad Cell Type Proportions near Megakaryocytes",
    x = "Distance (µm)", y = "Cell Proportion (%)", fill = "Broad class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


