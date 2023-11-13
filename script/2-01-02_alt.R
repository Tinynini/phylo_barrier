#library(tidyverse)

##########################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                   #
# Utilite == nettoyer la matrice et extraire les especes #
# Input == t_full_matrix.tsv                             #
# Output == sliced_all_species_clust.tsv et species.tsv  #
##########################################################

#### Ouverture de matrix.tsv & recuperation des donnees dans des dataframes ####
matrix <- read.csv('W:/Methode_Generic/output/matrix.tsv', header = TRUE, sep = ",")

rownames(matrix) <-(matrix[, 1])
matrix <- matrix[, -1]
matrix <- t(matrix)

#### Suppression des fichiers vides ####
keep_row <- rep(FALSE, nrow(matrix))

for (i in 1:nrow(matrix))
{
  keep_row[i] <- sum(matrix[i,]) > 0
}

matrix <- matrix[keep_row,]

#### Extraction et pretraitement des noms d'especes ####
species <- as.data.frame(sort(rownames(matrix)))
colnames(species) <- 'species'

matrix <- as.data.frame(matrix)

matrix %>%
  mutate(species = rownames(matrix), .before = matrix[, 1]) %>%
  arrange(species) %>%
  identity() -> matrix

err <- startsWith(species[, 'species'], 'UNV')

species[err, 'species'] <- str_replace(species[err, 'species'], pattern = "(.*)_(.*)_(.*)", replacement = "\\3\\_\\1\\_\\2")

matrix[, 1] <- species[, 1]

species %>%
  arrange(species) %>%
  identity -> species

matrix %>%
  arrange(species) %>%
  identity() -> matrix

#### Enregistrement de la dataframe dans le fichier sliced_all_species_clust.tsv ####
write.table(matrix, "W:/Methode_Generic/output/Dataframe/sliced_all_species_clust.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(species, "W:/Methode_Generic/output/Dataframe/species.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)