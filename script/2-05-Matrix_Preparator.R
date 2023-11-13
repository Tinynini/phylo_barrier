#library(tidyverse)

##############################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                       #
# Utilite == generer la 1ere matrice d absence/presence      #
# GenexNiveau et la liste des genes en vue du script suivant #
# Input == sliced_all_species_taxo.tsv                       #
# Output == Sliced_Matrix_species.tsv et uni_gene.tsv        #
##############################################################

#### Ouverture de sliced_all_species_taxo.tsv & recuperation des donnees dans une dataframe ####
all_species <- read_tsv('W:/Methode_Generic/output/Dataframe/sliced_all_species_taxo.tsv', show_col_types = FALSE) %>%
  as.data.frame()

keep_col <- rep(FALSE, ncol(all_species))
keep_col[1] <- TRUE

for (i in 2:ncol(all_species))
{
  keep_col[i] <- sum(all_species[, i]) > 0
}

all_species <- all_species[, keep_col]

#### Obtention & enregistrement de la matrice binaire GenexEspece ####

curr_matrix <- as.matrix(all_species[, -1])
rownames(curr_matrix) <- c(all_species[, 1])

write.table(curr_matrix, "W:/Methode_Generic/output/Matrice/Sliced_Matrix_species.tsv", sep = '\t', row.names = TRUE, col.names = TRUE)

#### Obtention & enregistrement de la liste des genes ####
uni_gene <- sort(colnames(all_species[, -1])) # On extrait les genes

uni_gene <- as.data.frame(uni_gene)
colnames(uni_gene) <- "qseqid"

write.table(uni_gene, "W:/Methode_Generic/output/Dataframe/uni_gene.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)