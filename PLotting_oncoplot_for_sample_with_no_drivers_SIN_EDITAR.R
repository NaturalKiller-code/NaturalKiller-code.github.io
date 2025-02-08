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

# Graficamos
oncoplot(maf = maf_sin_driver, top = 10)
