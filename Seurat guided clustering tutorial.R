#THe tutorial to learn how to do single cell RNA seq https://satijalab.org/seurat/articles/pbmc3k_tutorial.html




#First, have to install the packages of libraries that you need to use. So you call on them with function "instll.packages()". 
#If you add lines later, you need to select the whole line and run it because R sees only the one that you leave your mouse on

#installing dplyr package and related ones just in case, to be able to manipulate data (also contains ggplot2)

install.packages("tidyverse") #from https://dplyr.tidyverse.org/#installation 

#installing Seurat toolkit for single cell RNA-seq data manipulation and quality control analysis from CRAN, developed by Satija Lab and Collaborators.Requires compilation.
install.packages("Seurat") #from https://satijalab.org/seurat/articles/install.html

install.packages("patchwork") #to combine graphs from ggplot - package for graph making

#Now trying demo for "Seurat - Guided Clustering Tutorial" as practice.

library(dplyr) #attaching dplyr package to the file
#browseVignettes(package = "dplyr") #additional explanation on the library

library(Seurat)
library(patchwork)

# Load the PBMC dataset (has to be unpacked raw data with full path included or it doesn't work)
pbmc.data <- Read10X(data.dir = "C:/Users/Ivana/Documents/IMIM/Bioinfo/R practice/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) -> pbmc
pbmc
#Returns: An object of class Seurat 
#13714 features across 2700 samples within 1 assay 
#Active assay: RNA (13714 features, 0 variable features)

#The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
#These represent the selection and filtration of cells based on quality control (QC) metrics, data normalization and scaling, and the detection of highly variable features.

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Metadata gets stored in the Seurat Object that was created in line 29

# Visualize QC metrics as a violin plot.
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filter the cells with counts less then 200 (false detection of probably empty droplets) 
#& those of more then 2500 (becuse it's probably 2 or more cells in the same count)
#& those that have more then 5% of mitochondrial reads, originating probably from dying cells that are going through mitochondria mediated cell death (labeled with MT in the line 39)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normlizing the resuts.

#By default, we employ a global-scaling normalization method "LogNormalize" 
#that normalizes the feature expression measurements for each cell by the total expression, 
#multiplies this by a scale factor (10,000 by default), and log-transforms the result.
#Normalized values are stored in "pbmc[["RNA"]]@data".

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#these are default values for NormalizeData function, and R performs default if nothing is specified so same thing would be
#if we just wrote this "pbmc <- NormalizeData(pbmc)"
pbmc.data #to see the normalized data

#Next features need to be selected as representatives of cell subtypes, and it's easier for analysis if they are sorted by subset of
#features that show high variation between cells 
#from https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #this returns 2000 features per dataset, 
#but I guess you could specify more in the "nfeatures" part

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels. repel means that on the graph labels don't overlap and are readable. 

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #labels come from the data itself?
plot1 + plot2

#nice package from CRAN for cute labels: "ggrepel" ; also "identify()" let's you maualy place them
#(https://stats.stackexchange.com/questions/16057/how-do-i-avoid-overlapping-labels-in-an-r-plot)


#Now that there are labels, next thing is to scale the data by linearly transforming it.
#There are genes that differ in level of expression across cells, but this analysis is not interested in it since it's hard to
#distingush it from technical errors - and those we need to remove with normalization (this whole step)
#but there's a paper that says how to do it better and differentiate cell states (from the same lab) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
#also have vignette about it https://satijalab.org/seurat/articles/sctransform_vignette.html

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Now to reduce linear dimensions with PCA analysis https://www.sartorius.com/en/knowledge/science-snippets/what-is-principal-component-analysis-pca-and-how-it-is-used-507186
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Examine and visualize PCA results in a few different ways. PC_N means certain cell? Or cell type? NO, both cells and features are included...
#this explains PCA in context of scRNA seq https://hbctraining.github.io/scRNA-seq/lessons/05_normalization_and_PCA.html
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #compares cells across pincipal components and features, but orders them
#so that extremely different cells are at opposite side of spectrum

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) #this is for 15 PCs at once

#Figure out how many PC are significant - dimensions
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15) #comparing the distribution of p-values for each PC with a uniform distribution (dashed line)

ElbowPlot(pbmc) #a ranking of principle components based on the percentage of variance explained by each one 

#the two plots should show difference between PC that significantly differ from each other 
#AND the ones that do not - cut off needs to be assessed by user! Because rare subsets might not be expressed in a lot of numbers

#CLUSTERING - to figure out graphically different cell states, phenotypes 
#they use a special technique to link the cells that are most similar and then cluster them together in 2D

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#the results need to be requested, like this
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5) #the first 5 cells and their IDs that enable the clustering
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 50)

#To show data in 2D, multiple dimensions need to be reduced to the level of a graph " to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space"
#commonly with non-linear dimensional reduction, either tSNE or UMAP methods

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
# this is to visualize the calculated graph of UMAP
DimPlot(pbmc, reduction = "umap", label = "TRUE", repel = "TRUE")

#to save the graph so it's not necessary to perform all the code before it needs to be converted to object
saveRDS(pbmc, file = "C:/Users/Ivana/Documents/IMIM/Bioinfo/R practice/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/pbmc_tutorial_UMAP.rds")

#FINALLY  to find a distinct features in clusters as biomarkers

#to see the all the markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25) #this is the calculation
head(cluster2.markers, n = 5) #this is visualization in console

# find all markers distinguishing cluster 5 from clusters 0 and 3 but with limiting cell numbers that are analyzed to save time
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25) 
#need to have a cluster's markers as objects to compare them!!!
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#test differential expression in several ways https://satijalab.org/seurat/articles/de_vignette.html
#example ROC test (on a sale of 0-1, from random to perfect marker):
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE) #saved as object

#Visualization of Marker Expression in clusters

#violin plot: shows expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

#paints a particular feature (marker) on UMAP plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

#"DoHeatmap()" generates an expression heatmap for given cells and features. 
#In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

RidgePlot(object = "pbmc", features = "top10$gene", assay = "NULL") #https://satijalab.org/seurat/reference/ridgeplot sounds good, doesn't work

#more functions https://satijalab.org/seurat/reference/cellscatter https://satijalab.org/seurat/reference/dotplot

#you can assign cell types based on markers, just need to replace the labels - or rather: IDs, names
new.cluster.ids <- c("Naive CD4 T", "CD14 Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A Mono", "NK", "DC", "Platelet") #has to be in the same row
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "C:/Users/Ivana/Documents/IMIM/Bioinfo/R practice/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/pbmc3k_final_identified_subsets.rds")
