#Loading library
library(dndscv)
#Reading input data
acral_data <- read.csv("filtered_unique_snp_strelka_indels.csv", header=TRUE, sep='\t')
#Running dndscv
dndsout = dndscv(acral_data, refdb = "RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL)
#Showing significant genes
sel_cv = dndsout$sel_cv
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

##################################################################

library(maftools)
library(RColorBrewer)
#Reading mafs containing SNVs (folder contains only one maf per patient)
mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/SNP", pattern = "*.\\.csv$", full.names =TRUE)
#Merging  SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading mafs containing indels (folder contains only one maf per patient)
indel_mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/INDELS", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indek mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))
#Reading file containing clinical data
clindata = read.csv("Supplementary_Table_1.csv")


sample_summary_a <- getSampleSummary(indel)
sample_summary_b <- getSampleSummary(all_samples)
sample_summary_final <- getSampleSummary(all_alt)

# Identifica diferencias
unique_samples_a <- sample_summary_a$Tumor_Sample_Barcode
unique_samples_b <- sample_summary_b$Tumor_Sample_Barcode
unique_samples_final <- sample_summary_final$Tumor_Sample_Barcode

# Encuentra la muestra "extra"
extra_sample <- setdiff(unique_samples_final, union(unique_samples_a, unique_samples_b))
print(extra_sample)


###############
# Extraer los nombres de las muestras de cada archivo MAF (antes de fusionar)
muestras_a <- indel@data$Tumor_Sample_Barcode
muestras_b <- all_samples@data$Tumor_Sample_Barcode

# Verificar si hay nombres duplicados en cada archivo (antes de fusionar)
duplicados_a <- muestras_a[duplicated(muestras_a)]
duplicados_b <- muestras_b[duplicated(muestras_b)]

# Imprimir los duplicados si existen
print(duplicados_a)
print(duplicados_b)

# Verificar si hay muestras comunes entre los dos archivos (pueden ser diferentes por un sufijo, espacio, etc.)
comunes <- intersect(muestras_a, muestras_b)
print(comunes)

# Verificar las muestras que están en A pero no en B, y viceversa
solo_a <- setdiff(muestras_a, muestras_b)
solo_b <- setdiff(muestras_b, muestras_a)

# Imprimir las muestras únicas de cada archivo
print(solo_a)
print(solo_b)





###############

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Accent')
#Matching colors with categories
names(clincolors) = c("FEET", "HAND", "SUBUNGUAL")
clincolors = list(Site = clincolors)

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Site"), bgCol="white")


#Plotting other annotations

#Gender

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
#Matching colors with categories
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"), bgCol="white")


#Age

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"), bgCol="white")

#Ulceration

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
#Matching colors with categories
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Ulceration"), bgCol="white")

#Tumor_type
#For this annotation, tumor type was simplified in four categories: primary, metastasis, recurrence and lesion in transit
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"), bgCol="white")


#Stage

#Reading stage data as factors
clindata$Stage <- as.factor(clindata$Stage)
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("NRAS","BRAF","KIT","NF1","HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","CSN2","PDGFRA","MED12","ZNF793","MDM1","KRTAP4-16","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")





##########################

## Visualizing mutations in specific genes
#Lollipop plots
lollipopPlot(maf=all_alt, gene="NRAS",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="BRAF",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="NF1",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
lollipopPlot(maf=all_alt, gene="KIT",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

#############


## PLotting oncoplot for sample with no drivers

#Loading packages
library(maftools)
library(RColorBrewer)
#Reading MAFs with SNV data
mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/SNP", pattern = "*.\\.csv$", full.names =TRUE)
#Merging  SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading MAFs with indel data
indel_mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/INDELS", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indel mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))
#Reading file containing clinical data
clindata = read.csv("Supplementary_table_1.csv")

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Accent')
names(clincolors) = c("FEET", "HAND", "SUBUNGUAL")
clincolors = list(Site = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Site"), bgCol="white")


#Gender

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=TRUE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"), bgCol="white", SampleNamefontSize=0.8)


#Age

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"), bgCol="white", SampleNamefontSize=0.8)


#Ulceration

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Ulceration"), bgCol="white")


#Tumor_type

#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"), bgCol="white")


#Stage

#Choosing colors
clindata$Stage <- as.factor(clindata$Stage)
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")




#################
#Reading MAFs with SNV data
mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/SNP", pattern = "*.\\.csv$", full.names =TRUE)
#Merging  SNV mafs
all_samples = merge_mafs(mafs = mafs)
#Reading MAFs with indel data
indel_mafs = list.files(path = "/Users/emilianoferrorodriguez/Desktop/LCGEJ/Cancer/Muestras/INDELS", pattern = "*.\\.maf$", full.names =TRUE)
#Merging indel mafs
indel = merge_mafs(mafs = indel_mafs)
#Merging SNV and indel mafs
all_alt = merge_mafs(mafs= c(indel, all_samples))


#Plotting lollipo plots

#HRAS

lollipopPlot(maf=all_alt, gene="HRAS",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

#SPHKAP

lollipopPlot(maf=all_alt, gene="SPHKAP",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)

#POU3F3

lollipopPlot(maf=all_alt, gene="POU3F3",  AACol = 'HGVSp_Short', labelPos = 'all', showMutationRate=TRUE)
