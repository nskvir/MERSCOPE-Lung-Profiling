#updated Giotto Lung Analysis Script


#pak::pkg_install("drieslab/Giotto")

library(cli)
library(pak)
library(Giotto)
library(bit64)
library(arrow)
library(future.apply)


root_dir <- "/path/to/data/directory/"

runs <- list(
    A = file.path(root_dir, "202403242054_VA00340-BUMC-CellBoundary-GroupA_VMSC01801", "region_0"),
    B = file.path(root_dir, "202403242055_VA00340-BUMC-CellBoundary-GroupB_VMSC03501", "region_0")
)


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

#create giotto object for a sample
g <- giotto()
g <- setGiotto(g, polygons)
g <- setGiotto(g, gpoints)
g <- setGiotto(g, imgs)


###########################
# This portion of the script is just for tesselating the image to manually look for megs

# Add hex data to object - can shrink hexagons to shape_size = 50 as needed to zoom in
hex <- tessellate(extent = ext(imgs$PolyT$z0), shape = "hex", shape_size = 100, name = "hex100")
g <- setGiotto(g, hex)

#calculate hex overlap information
g <- calculateOverlap(g, spatial_info = "hex100", feat_info = "rna")
g <- overlapToMatrix(g, poly_info = "hex100", feat_info = "rna")
g <- filterGiotto(g, spat_unit = "hex100", feat_det_in_min_cells = 1,
                  min_det_feats_per_cell = 1, expression_threshold = 1)
g <- normalizeGiotto(g, spat_unit = "hex100")

# signature metafeature ####
# create a summary score for all the genes in the mk signature so that
# it is easier to visualize where the values are concentrated

# setup mk gene signature
mk <- c(
    "Gp1ba", "Itga2b", "Mpl", "Pf4", "Tubb1", "Vwf"
)
anno <- rep("mk_sig", 6)
names(anno) <- mk
# rescale each gene to range between 0 and 1
# find the sum value of the rescaled genes per cell as the overall
# metafeature score to report for that cell
g <- createMetafeats(g, spat_unit = "hex100", name = "rescaled_sum_sig",
                     feat_clusters = anno, rescale_to = c(0, 1),
                     stat = "sum")
# m
g <- createMetafeats(g, spat_unit = "hex100", name = "scaled_mean_sig",
                     feat_clusters = anno,
                     expression_values = "scaled",
                     stat = "mean")

# increase resampling size so images are clearer
options("giotto.plot_img_max_sample" = 1e7)



### Plots used in Figure 3A

spatInSituPlotPoints(
    g, show_image = T, image_name = "DAPI_z5", spat_unit = "hex100",
    spat_enr_names = "rescaled_sum_sig", polygon_fill = "mk_sig",
    polygon_fill_gradient_style = "sequential",
    feats = list(rna = mk), point_size = 0.2,
    save_plot=F,
    save_param = list(save_name = file.path('Figures/GroupA_Hex_rescaled_sum'),
    dpi = 600)
)

spatInSituPlotPoints(
    g, show_image = T, image_name = "DAPI_z5", spat_unit = "hex100",
    spat_enr_names = "scaled_mean_sig", polygon_fill = "mk_sig",
    polygon_fill_gradient_style = "divergent",
    polygon_fill_gradient = c("cyan", "blue", "black", "orange", "yellow"),
    polygon_fill_gradient_midpoint = 0,
    feats = list(rna = mk), point_size = 0.2,
    save_plot=F,
    save_param = list(save_name = file.path('Figures/GroupA_Hex_scaled_mean'),
    dpi = 600)

)

# Zoom in on potential mk areas ####

# Test a specified subset/FOV using x,y coordinates
mini <- subsetGiottoLocs(
    g, spat_unit = ":all:",
    x_min = 4000, x_max = 5000, y_min = 4000, y_max = 5000, z_min = NULL, z_max = NULL,
    poly_info = list(hex100 = "hex100")
    )


# Plot DAPI (or Cellbound1/2/3, PolyT, etc) images
spatInSituPlotPoints(
    mini, show_image = T, image_name = "DAPI_z4",
    feat_type = "rna",
    spat_unit = "hex100",
    show_polygon = F,
    use_overlap = F,
    point_size = 0.5,
    save_plot=F,
    save_param = list(save_name = file.path('Figures/GroupA_FOV1_DAPI_z4'),
                      dpi = 600)
)

### Plot images in conjunction with transcripts (utilized for Figures 3B and 3C)
spatInSituPlotPoints(
    mini, show_image = T, image_name = "Cellbound1_z3", 
    show_polygon = F,
    use_overlap = F,
    spat_unit = "hex100",
    spat_enr_names = "rescaled_sum_sig",
    polygon_fill = "mk_sig",
    polygon_alpha = 0.2,
    polygon_fill_gradient_style = "sequential",
    feats = list(rna = mk),
    point_size = 0.5,
    save_plot=F,
    save_param = list(save_name = file.path('Figures/GroupA_FOV1_x123_456_y123_456_Cellbound1_z3'),
                      dpi = 600)
)
############################

