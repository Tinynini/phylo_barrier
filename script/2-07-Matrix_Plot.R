#library(tidyverse)

###################################################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                                            #
# Utilite == generer les barplot de presence d un gene donne avec les 15 matrices pseudo_binaires #
# Input == taxo_species.tsv, uni_gene.tsv et les 15 matrices pseudo-binaires                      #
# Output == les 15 barplots (en FR et EN)                                                         #
###################################################################################################

#### Ouverture de taxo_species.tsv et uni_gene.tsv & recuperation des donnees ####

species <- read_tsv('W:/Methode_Generic/output/Dataframe/taxo_species.tsv', col_types = 'cccccc') %>% 
  as.data.frame()

uni_gene <- read_tsv('W:/Methode_Generic/output/Dataframe/uni_gene.tsv', col_types = 'c') %>% 
  as.data.frame()

#### Choosing one particular gene by it's index
gene_index = 10

#### Obtention de la liste des fichiers contenant nos matrices binaires et pseudo-binaires ####
# On se rend dans le bon emplacement puis on recherche la parterne de nom de fichier 'Sliced_Matrix_.*.tsv' qui est commune aux 2 types de matrice 
all_matrix <- list.files(path = 'W:/Methode_Generic/output/Matrice', pattern = 'Sliced_Matrix_.*.tsv', full.names = TRUE) 
n_matrix <- length(all_matrix) # On recupere le nombre de fichier contenus dans notre liste de fichier 

matrix_name <- str_replace(all_matrix, '(.*)(Matrix)_(.*).(tsv)', '\\3') # On recupere les noms de matrices a partir des noms de fichiers 

for (i in 1:n_matrix) # Permet de parcourir les matrices une par une
{
  if (grepl('(.*)_(.*)', matrix_name[i]) == TRUE) # On ne fait ce plot que pour les matrices pseudo-binaires
  { 
    gene_matrix <- read.csv(file = all_matrix[i], header = TRUE, sep = "\t") # On ouvre la matrice depuis la liste de fichier
    gene_matrix <- as.matrix(gene_matrix)
    
    #### Barplot de la presence d un gene donne ####
    # Pour definir les noms et destinations de fichiers pour l enregistrement
    deb_fr <- "W:/Methode_Generic/output/Plot/Matrice_plot/Barplot_Pres_Gene/rep/FR/Partage_inter-"
    fin_fr <- "_fr.png"
    deb_en <- "W:/Methode_Generic/output/Plot/Matrice_plot/Barplot_Pres_Gene/rep/EN/Partage_inter-"
    fin_en <- "_en.png"
    
    # Pour definir les titres de barplots
    debut <- "Partage inter-"
    fin <- " de "
    start <- " inter-"
    end <- " sharing"
    
    # N.B. : Il suffit de changer l index dans gene_matrix et uni_gene dans les titres de graphs pour tester un autre gene
    to_set <- which(gene_matrix[, gene_index] != 0) # On isole les lignes pour lesquelles le gene matche
    m <- gene_matrix[c(to_set), gene_index] # On extrait lesdites lignes
    
    png(str_glue("{deb_fr}{matrix_name[i]}{fin_fr}"), height = 1017, width = 1920, pointsize = 20)
    barplot(m, main = str_glue("{debut}{matrix_name[i]}{fin}{uni_gene[gene_index,]}")) # Barblot de la presence du gene donne dans la matrice
    dev.off()
    png(str_glue("{deb_en}{matrix_name[i]}{fin_en}"), height = 1017, width = 1920, pointsize = 20)
    barplot(m, main = str_glue("{uni_gene[gene_index,]}{start}{matrix_name[i]}{end}")) # Barblot de la presence du gene donne dans la matrice
    dev.off()
  }
}