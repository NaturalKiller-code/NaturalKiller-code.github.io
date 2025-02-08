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

# Listamos los genes Driverss
genes_drivers <- c("NRAS", "BRAF", "NF1", "KIT")

# Quitamos los genes drivers
gene_Summary <- getGeneSummary(all_alt)
genes_mutados <- gene_Summary$Hugo_Symbol

# Pedimos todas las muestras
muestras_completas <- getSampleSummary(all_alt)
muestras_completas$Tumor_Sample_Barcode

# Revisar los datos clínicos en cada archivo
head(all_samples@clinical.data)
head(indel@clinical.data)



## Buscamos las muestras sin drivers

# Creamos lista de genes drivers
genes_drivers <- c("NRAS", "BRAF", "NF1", "KIT")

# Filtramos las muestras que contengan las mutaciones en la lista gene_drivers
mutaciones_driver <- subsetMaf(all_alt, genes = genes_drivers)

# Pedimos todas las muestras
muestras_completas <- getSampleSummary(all_alt)$Tumor_Sample_Barcode

# Pedimos todas las muestras sin genes drivers
muestras_sin_drivers <- setdiff(muestras_completas, mutaciones_driver@clinical.data$Tumor_Sample_Barcode)
print(muestras_sin_drivers)


## Creamos los Oncoplots
# Creamos MAF con las muestras sin drivers
maf_sin_driver <- subsetMaf(all_alt, tsb = muestras_sin_drivers)

# Creamos un archivo csv para la muestra 
# Acceder a la tabla de mutaciones asociada a las muestras sin mutaciones en genes drivers
tabla_mutaciones <- maf_sin_driver@data

# Cargamos los datos en un archivo CSV
write.csv(tabla_mutaciones, "muestras_sin_drivers.csv", row.names = FALSE)

# Verificamos el archivo csv
datos_muestras_sin_drivers <- read.csv("muestras_sin_drivers.csv")


# Graficamos
oncoplot(maf = maf_sin_driver, top = 10)

#Gender

# With other genes, without drivers
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=TRUE, clinicalFeatures= c("Gender"), bgCol="white", SampleNamefontSize=0.8)


# Ignore genes
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 3,name = 'Dark2')
names(clincolors) = c("female", "male", "NA")
clincolors = list(Gender = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=TRUE, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Gender"), bgCol="white", SampleNamefontSize=0.8)




#Age

#Plotting oncoplot
oncoplot(all_alt, draw_titv = TRUE, annotationDat=clindata, showTumorSampleBarcodes=FALSE, genes = c("HRAS","POU3F3","SPHKAP","RDH5","DLX6","SPRED1","CDKN3","SMCO3","PDGFRA","MED12","ZNF793","MDM1","PAXIP1","TP53","KRAS"), genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Age"), bgCol="white", SampleNamefontSize=0.8)


#Ulceration

# With other genes, without drivers
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, clinicalFeatures= c("Ulceration"), bgCol="white")

# Ignore genes 
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Set1')
names(clincolors) = c("yes", "no", "")
clincolors = list(Ulceration = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE,genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2")
, clinicalFeatures= c("Ulceration"), bgCol="white")





#Tumor_type

# With other genes without drivers
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, clinicalFeatures= c("Tumor_type"), bgCol="white")


# Ignore genes
#Choosing colors
clincolors = RColorBrewer::brewer.pal(n = 4,name = 'Paired')
names(clincolors) = c("primary", "metastasis", "Recurrence","Lesion_in_transit")
clincolors = list(Tumor_type = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Tumor_type"), bgCol="white")


#Stage


# Ignore Genes
#Choosing colors
clindata$Stage <- as.factor(clindata$Stage)
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, genesToIgnore = c("TTN","OBSCN", "MUC16", "USH2A","SYNE2"), clinicalFeatures= c("Stage"), bgCol="white")



# With another genes without drivers
clindata$Stage <- as.factor(clindata$Stage)
clincolors = RColorBrewer::brewer.pal(n = 5,name = 'BuPu')
names(clincolors) = c("0", "1", "2","3","4")
clincolors = list(Stage = clincolors)
#Plotting oncoplot
oncoplot(maf = maf_sin_driver, draw_titv = TRUE, annotationDat=clindata, annotationColor = clincolors, showTumorSampleBarcodes=FALSE, clinicalFeatures= c("Stage"), bgCol="white")





##########################
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






#########################

#Loading library
library(dndscv)
#Reading input data
acral_data <- read.csv("filtered_unique_snp_strelka_indels.csv", header=TRUE, sep='\t')

# Pedimos todas las muestras sin genes drivers
muestras_sin_drivers <- setdiff(muestras_completas, mutaciones_driver@clinical.data$Tumor_Sample_Barcode)
print(muestras_sin_drivers)
 
snp_sin_drivers <- acral_data[acral_data$sampleID %in% datos_muestras_sin_drivers$Tumor_Sample]

#Running dndscv
dndsout = dndscv(acral_data, refdb = "RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=NULL)
#Showing significant genes
sel_cv = dndsout$sel_cv
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)






