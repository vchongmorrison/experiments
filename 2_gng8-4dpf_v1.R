### Adapted from https://github.com/hbctraining/scRNA-seq/tree/master/lessons
## Load Seurat and other required packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(AnnotationHub) # needs "ensembldb" installed

##### Create Seurat object from Cell Ranger output #####

## Import the raw, unfiltered matrix
# How to read in 10X data for a single sample (output is a sparse matrix)
counts <- Read10X(data.dir = "./raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
gng8 <- CreateSeuratObject(counts = counts,
                           min.features = 100) # only consider 'cells' with > 100 genes detected

# Explore the metadata
View(gng8@meta.data)

# Add number of genes per UMI for each cell to metadata
gng8$log10GenesPerUMI <- log10(gng8$nFeature_RNA) / log10(gng8$nCount_RNA)

# Compute proportion of tx mapping to mitochondrial genes

gng8$mitoRatio <- PercentageFeatureSet(object = gng8, pattern = "^mt-")
gng8$mitoRatio <- gng8@meta.data$mitoRatio / 100 # compute ratio instead of percentage

# Create metadata dataframe for further messing about without affecting the original metadata.
metadata <- gng8@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                       nUMI = nCount_RNA,
                                       nGene = nFeature_RNA)

# Create sample column
metadata$sample <- "gng8-4dpf"

# Add metadata back to Seurat object
gng8@meta.data <- metadata
View(gng8@meta.data)

# Create .RData object to load at any time
save(gng8, file="./gng8_filtered_100.RData")

##### Quality Check of data #####

# Cell counts
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells with min=100 genes")

# UMI counts per cell
# Visualize the number UMIs/transcripts per cell (should generally be above 500)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500)

# Genes detected per cell
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# UMIs vs genes detected
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Mitochondrial counts ratio
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Complexity
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

##### Filtering for good quality cells #####

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_gng8 <- subset(x = gng8, 
                        subset= (nUMI >= 180) & 
                          (nGene >= 150) & 
                          (log10GenesPerUMI > 0.80) & 
                          (mitoRatio < 0.20)) 

nrow(filtered_gng8@meta.data) # 1224 cells remain

##### Filtering genes that are only expressed in X or more cells #####

# Extract counts
counts <- GetAssayData(object = filtered_gng8, slot = "counts") 
nrow(counts) # 25108 genes

# Output a logical vector for every gene on whether there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 1 # going lenient and accept if expressed in ONE cell 

# Only keeping those genes expressed in more than X cells
filtered_counts <- counts[keep_genes, ] 
nrow(filtered_counts) # 20362 genes remain

# Reassign to filtered Seurat object
filtered_gng8 <- CreateSeuratObject(filtered_counts, meta.data = filtered_gng8@meta.data) 
nrow(filtered_gng8@meta.data) # still 1224 cells - OK

# Save filtered subset to new metadata to re-visualise QC metrics as above
metadata_clean <- filtered_gng8@meta.data

# Create .RData object to load at any time
save(filtered_gng8, file="./gng8_filtered.RData")

##### Check cell x gene matrix is in order #####
# Load gng8_filtered.RData
counts <- GetAssayData(object = filtered_gng8, slot = "counts")
nrow(counts) # 20362 genes
ncol(counts) # 1224 cells

##### Obtain cell cycle genes #####

# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:

cc_file <- RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Danio_rerio.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Danio rerio", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Save annotations file for later use
write.csv(annotations, file='annotations.csv')

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

###### Score cells for cell cycle #####

# Roughly normalize the counts (not as accurate as SCtransform which will do later for Clustering)
seurat_phase <- NormalizeData(filtered_gng8)

seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)            

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase") # can't see any obvious differences due to cell cycle phase - don't regress

##### SCTransform (normalise, scale, find variable genes) - single sample #####
filtered_gng8 <- SCTransform(filtered_gng8, vars.to.regress = "mitoRatio", verbose = FALSE)
filtered_gng8 <- CellCycleScoring(filtered_gng8, g2m.features = g2m_genes, s.features = s_genes) # add this into object

# Save as RData file to load anytime for subsequent clustering and visualisation.
save(filtered_gng8, file="./gng8_filtered_SCT.RData")

##### Clustering cells based on top PCs (metagenes) #####

## Run PCA
filtered_gng8 <- RunPCA(object = filtered_gng8) # 50 PCs by default

## Identify significant PCs
# Explore heatmap of PCs
DimHeatmap(filtered_gng8, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = filtered_gng8[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = filtered_gng8, 
          ndims = 40) # Inflexion around 8-10ish PC?

## Cluster the cells
# Determine the K-nearest neighbor graph
filtered_gng8 <- FindNeighbors(object = filtered_gng8, 
                               dims = 1:40) # first 40 dimensions

# Determine the clusters for various resolutions                                
filtered_gng8 <- FindClusters(object = filtered_gng8,
                              resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)) # recommended "granularity" for 3-5K cells

# Explore resolutions
filtered_gng8@meta.data %>% 
  View()

# Save as RData file to load anytime for subsequent DE analysis.
save(filtered_gng8, file="./gng8_filtered_SCT_clusters.RData")

# Assign identity of clusters for further analysis
Idents(object = filtered_gng8) <- "SCT_snn_res.0.4" # using 0.4 resolution

# Plot the UMAP (or tsne, pca)

# Run UMAP
filtered_gng8 <- RunUMAP(filtered_gng8, # or RunTSNE
                         dims = 1:40,
                         reduction = "pca")

DimPlot(filtered_gng8,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

##### Exploration of PC metrics #####

## Segregation of clusters by sample
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(filtered_gng8, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
write.table(n_cells, file='./cluster04_ncells.txt') # save this information

## Segregation of clusters by cell cycle phase
# Explore whether clusters segregate by cell cycle phase
DimPlot(filtered_gng8,
        label = TRUE, 
        split.by = "Phase")  + NoLegend() # nothing jumps out...

## Segregation of clusters by various sources of uninteresting variation
# Determine metrics to plot present in filtered_gng8@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(filtered_gng8, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

## Exploration of PCs driving different clusters
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16), # Look at top 16 PCs
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(filtered_gng8, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(filtered_gng8, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results 
print(filtered_gng8[["pca"]], dims = 1:5, nfeatures = 20) # 20 features of top 5 PCs

##### Exploring known cell type markers #####
DimPlot(object = filtered_gng8,
        reduction = "umap", 
        label = TRUE) + NoLegend()

# Select the raw RNA counts slot to be the default assay
DefaultAssay(filtered_gng8) <- "RNA"

# Normalise RNA data for visualization purposes
filtered_gng8 <- NormalizeData(filtered_gng8, verbose = FALSE)

# Plot known markers - some examples
FeaturePlot(filtered_gng8, 
            reduction = "umap", 
            features = c("pvalb5", "calb2a", "epcam", "s100z"), # olfactory
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(filtered_gng8, 
            reduction = "umap", 
            features = c("kctd12.1", "nrp1a", "kctd8", "slc18a3b"), # dHb_left, dHb_right
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

##### Identification of all markers for each cluster #####
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = filtered_gng8, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

# Add gene annotations
annotations <- read.csv('./annotations.csv', header=TRUE, row.names=1)

# Combine markers with gene descriptions (single sample)
markers_ann <- markers %>% left_join(y = unique(annotations[, c("gene_name", "description")]),
                                     by = c("gene" = "gene_name"))

View(markers_ann)
write.csv(markers_ann, file='./markers.csv')

##### Evaluating marker genes #####
markers_ann <- read.csv('./markers.csv', header=TRUE, row.names=1)

# Extract top 10 markers per cluster (if single sample)
top10 <- markers_ann %>% 
  group_by(cluster) %>% 
  top_n(n = 10, 
        wt = avg_log2FC)

# Visualize top 20 markers per cluster
View(top10)

# Plot heatmap of top X markers per cluster
# Select the scaled counts slot to be the default assay
DefaultAssay(filtered_gng8) <- "SCT" # use scaled (counts) data, *n.b. only has variable genes

heatmap_features <- top10$gene
DoHeatmap(filtered_gng8, assay = "SCT", features = heatmap_features, size = 3)

# If want a 'targeted' heatmap...
heatmap_features <- c("cachd1", "lrp5", "lrp6", "fzd7a", "fzd7b")
cells0 <- WhichCells(filtered_gng8, idents="0") # cluster of interest
cells1 <- WhichCells(filtered_gng8, idents="1") # cluster of interest
cells3 <- WhichCells(filtered_gng8, idents="3") # cluster of interest
cells <- c(cells0,cells1,cells3)
DoHeatmap(filtered_gng8, assay = "SCT", features = heatmap_features, size = 3, cells=cells)

##### Visualisations of marker genes expression #####
# Load gng8_filtered_SCT_clusters.RData if required.
# Assign identity of clusters
Idents(object = filtered_gng8) <- "SCT_snn_res.0.4" # using 0.4 resolution

DimPlot(object = filtered_gng8, 
        reduction = "umap", 
        label = TRUE) + NoLegend() # check 

# Select the raw RNA counts slot to be the default assay
DefaultAssay(filtered_gng8) <- "RNA"

# Normalize RNA data for visualization purposes
filtered_gng8 <- NormalizeData(filtered_gng8, verbose = FALSE)

# Ridge plots - from ggridges. Visualise single cell expression distributions in each cluster.
filtered_gng8$sample <- sample(c("gng8-4dpf"), size = ncol(filtered_gng8), replace = TRUE)
features <- c("kctd12.1", "kctd12.2", "nrp1a", "kctd8", "slc18a3b", "pou4f1", "kiss1", "aoc1", "gfap", "cxcr4b", "dbx1b") 

RidgePlot(filtered_gng8, features = features, ncol = 6)

# Violin plot
VlnPlot(object = filtered_gng8, 
        features = c("kctd12.1", "kctd12.2", "kctd8", "slc18a3b", "cxcr4b", "cachd1"))
VlnPlot(object = filtered_gng8, 
        features = c("kctd12.1", "kctd12.2", "kctd8", "slc18a3b", "cachd1", "lrp5", "lrp6", "fzd7a", "fzd7b"))
VlnPlot(object = filtered_gng8, 
        features = c("gng8", "elavl3"))

# Visualise co-expression of two features simultaneously
FeaturePlot(filtered_gng8, features = c("cachd1", "lrp5"), 
            blend = TRUE, cols = c('#1b9e77','#d95f02','#7570b3'), # colorblind-safe colorbrewer scale
            order=TRUE,
            pt.size=2,
            label=TRUE,
            repel=TRUE,
            reduction="umap") 

VlnPlot(filtered_gng8, features = "mitoRatio") # check status of cells in each cluster
VlnPlot(filtered_gng8, features = "nUMI")
VlnPlot(filtered_gng8, features = "log10GenesPerUMI")

##### Rename all identities ####
filtered_gng8 <- RenameIdents(object = filtered_gng8, 
                              "0" = "0_kctd12.1/2_kctd8_slc18a3b_pou4f1_cachd1?",
                              "1" = "1_cachd1? fzd7a/b? lrp5?",
                              "2" = "2_Olf1",
                              "3" = "3_Mystery Cluster",
                              "4" = "4_Olf2",
                              "5" = "5_Olf3",
                              "6" = "6_Olf4",
                              "7" = "7_kiss1_pou4f1_slc18a3b_kctd8_aoc1_kctd12.2",
                              "8" = "8_kiss1_pou4f1_slc18a3b_kctd8",
                              "9" = "9_cxcr4b_highTx_lowComplexity",
                              "10" = "10_cxcr4b_lowTx_highComplexity")


# (Re-)Plot the UMAP for presentations etc.
baseplot <- DimPlot(object = filtered_gng8, 
                    reduction = "umap", 
                    label = FALSE,
                    label.size = 3,
                    repel = TRUE)
baseplot+DarkTheme()

# Save final R object
write_rds(filtered_gng8,
          file = "./filtered_gng8_labelled.rds")       