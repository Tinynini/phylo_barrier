#library(tidyverse)

####################################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                             #
# Utilite == recuperer la taxonomie d un maximum d espece via des joins successifs #
# Input == sliced_all_species_clust.tsv, species_tsv et Parsed_taxonomy.tsv (V1)   #
# Output == sliced_all_species_clust.tsv et taxo_species.tsv                       #
####################################################################################

#### Ouverture de Parsed_taxonomy.tsv et de sliced_all_species_clust.tsv & recuperation des donnees ####
Parsed_taxonomy <- read_tsv('W:/Methode_Generic/output/Table_taxonomie/Parsed_taxonomy.tsv', col_types = "ccccccc")
Parsed_taxonomy <- Parsed_taxonomy[-c(11930, 11932, 10249, 14172, 16358, 16359),] # Suppression preventive de certaine lignes de la table de taxonomie pour eviter l apparition de certains doublons 

all_species <- read_tsv('W:/Methode_Generic/output/Dataframe/sliced_all_species_clust.tsv', show_col_types = FALSE) %>%
  as.data.frame()

species <- read_tsv('W:/Methode_Generic/output/Dataframe/species.tsv', show_col_types = FALSE) %>%
  as.data.frame()

#### Fonction servant a effectuer un traitement supplementaire pour les noms d especes qui matchent 'pattern_1' ####
special_treat <- function(df, j, pattern_1, pattern_2, replacement)
{
  if (grepl(pattern_1, df[j,'species']) == TRUE)
  {
    df[j,'species'] <- str_replace(df[j,'species'], pattern_2, replacement)
  }
  
  return(df)
}

#### Fonction servant a supprimer de maniere systematique certain types de doublons post-joint ####
Genus_cleaner <- function(df, level, doublon, ref)
{
  double1 <- which(df[, level] %in% doublon)
  double2 <- which(df[c(double1), 'Family'] != ref)
  gene_double <- df[c(double1),]
  df[c(double1),] <- gene_double[-c(double2),]
  
  df <- unique(df)
}

#### Fonction (utilisant celle ci-avant) servant a gerer l ensemble des doublon crees lors du join au niveau des genus ####
Genus_cleaning <- function(df)
{ # /!\ Les cas a traiter dependent des especes traitees et de la table de taxonomie utilisee comme ref !
  df <- Genus_cleaner(df, 'species', c('Clostridium aldenense', 'Clostridium clostridioforme'), 'Lachnospiraceae')
  df <- Genus_cleaner(df, 'species', 'Clostridium difficile', 'Peptostreptococcaceae')
  df <- Genus_cleaner(df, 'species', 'Clostridium phoceensis', 'Acutalibacteraceae')
  df <- Genus_cleaner(df, 'Genus', 'Ruminococcus', 'Lachnospiraceae')
  df <- Genus_cleaner(df, 'Genus', 'Eubacterium', 'Eubacteriaceae')
  df <- Genus_cleaner(df, 'Genus', 'Mycoplasma', 'Mycoplasmataceae')
  df <- Genus_cleaner(df, 'Genus', 'Ruminiclostridium', 'Acetivibrionaceae')
  df <- Genus_cleaner(df, 'Genus', 'Phormidium', 'Geitlerinemaceae')
  df <- Genus_cleaner(df, 'Genus', 'Spirochaeta', 'Spirochaetaceae')
  df <- Genus_cleaner(df, 'Genus', 'Bacteroides' , 'Bacteroidaceae')
  df <- Genus_cleaner(df, 'Genus', 'Desulfotomaculum' , 'Desulfotomaculaceae')
  df <- Genus_cleaner(df, 'Genus', 'Paenibacillus' , 'Paenibacillaceae')
  df <- Genus_cleaner(df, 'Genus', 'Rhodospirillum' , 'Rhodospirillaceae')
  
  double1 <- grep('Clostridium', df[, 'Genus'])
  gene_Clos <- df[c(double1),]
  double2 <- which(gene_Clos[, 'species'] %in% c('Clostridium aldenense', 'Clostridium clostridioforme', 'Clostridium difficile', 'Clostridium phoceensis'))
  gene_Target <- gene_Clos[-c(double2),]
  
  df <- Genus_cleaner(df, 'species', gene_Target[, 'species'], 'Clostridiaceae')
}

#### Fonction servant a recuperer manuellement la taxonomie des especes qui matchent 'wanted' ####
except_treat <- function(df, wanted, level_1, level_2, level_3, rep_1, rep_2, rep_3)
{ # On est oblige de traiter chaque niveau d interet separement autrement un decalage systematique a lieu a chaque nouvelle ligne traitee
  ex <- which(df[, 'species'] %in% wanted)
  df[ex, level_1] <- rep_1
  df[ex, level_2] <- rep_2
  df[ex, level_3] <- rep_3
  # On prevoit donc le traitement de 3 niveaux distincts (on n a pas besoin de plus heureusement)
  return(df)
} 

#### Preparation des donnees en vue d un 1er join au niveau des especes ####
# On modifie la nomenclature des noms d especes en vue du join 
species[, 'species'] <- str_replace(species[, 'species'], pattern = '(.*)_(.*)', replacement = "\\1\\ \\2")

for (j in 1:nrow(species)) # Certains noms d especes necessitent un traitement supplementaire
{ # /!\ Les cas a traiter dependent des especes traitees !
  species <- special_treat(species, j, 'ORG.', "(.*)_(.*) (.*)(.)", "\\1\\ \\2\\_\\3")
  species <- special_treat(species, j, 'CONTAM.', "(.*)_(.*) (.*)(.)", "\\1\\ \\2\\_\\3")
  species <- special_treat(species, j, 'symbiont', "(.*) (.*)", "\\2\\ \\1")
  species <- special_treat(species, j, 'Bacterium', "(.*) (.*)", "\\2\\ \\1")
}

#### 1er join au niveau des especes --> consequence : ajout de 6 nouvelles colonnes ('Genus' a 'Domain') ####

species <- left_join(species, Parsed_taxonomy, by = c('species' = 'Species'))

#### Preparation des donnees en vue des 3 joins successifs a venir ####
level_name <- unlist(colnames(Parsed_taxonomy[, c(1:4)]))
level_name[1] <- 'species'
# On liste les suffix des especes 'bacterium' pour pouvoir les traiter selon le niveau taxonomique a partir duquel elles sont referencees
suffix <- c('', '(.*) (bacterium)', '(.*)(ae) (bacterium)', '(.*)(ales) (bacterium)')
cond <- c('', FALSE, TRUE, TRUE) # On adapte la condition de traitement selon le niveau taxonomique
NA_level <- is.na(species[, level_name[2]]) # Extraction des especes qui n ont pas pu etre matchees lors du 1er join 
# /!\ Les suffix a traiter et les conditions associees dependent des especes traitees !
#### joins successifs au niveau des genus puis des familles et enfin des ordres ####
for (i in 2:4)
{
  # On supprime le niveau traite precedement de la table de taxonomie en vue du nouveau join
  Parsed_taxonomy <- Parsed_taxonomy[,-c(1)] 
  # Traitement specifiques des especes 'bacterium' qui ne peuvent etre matchees qu au niveau i 
  for (k in 1:nrow(species))
  {
    if (NA_level[k] == TRUE & grepl(suffix[i], species[k, 'species']) == cond[i])
    {
      species[k, level_name[i]] <- str_replace(species[k, 'species'], '(.*) (.*)', '\\1')
    }
  }
  # Creation d une nouvelle dataframe contenant uniquement les especes non matchee lors du join precedent
  gene_level <- as.data.frame(species[c(NA_level), c(1:i)])
  colnames(gene_level) <- colnames(species[, c(1:i)])
  gene_level <- left_join(gene_level, Parsed_taxonomy, by = NULL) # Join au niveau i
  gene_level <- unique(gene_level)
  # Traitement visant a supprimer les nombreux doublons generes par le join au niveau des genus ####
  if (i == 2)
  {
    gene_level <- Genus_cleaning(gene_level)
  }
  # Remplacement dans notre dataframe initiale des lignes associees aux especes non-matchees par celles de la nouvelle dataframe 
  less_NA_level <- which(species[, level_name[i - 1]] %in% gene_level[, level_name[i - 1]])
  species[c(less_NA_level),] <- gene_level
  NA_level <- is.na(species[, level_name[i]])
}

#### Traitement direct (hors join) de cas particuliers d especes (principalement des 'bactérium') trop isolees pour faire l objet d'un join ####
# /!\ Les cas a traiter dependent des especes traitees et de la table de taxonomie utilisee comme ref !
species <- except_treat(species, 'Bacillus bacterium', c('Genus', 'Family', 'Order'), 'Class', 'Phylum', c('Bacillus', 'Bacillaceae', 'Bacillales'), 'Bacilli', 'Firmicutes')
species <- except_treat(species, 'Lachnospiraceae oral', 'Genus', 'Family', 'Order', NA, 'Lachnospiraceae', 'Lachnospirales')
species <- except_treat(species, 'Sphingomonas.like bacterium', 'Genus', 'Family', 'Order', 'Sphingomonas.like', NA, NA)
species <- except_treat(species, 'Corynebacterium.like bacterium', 'Genus', 'Family', 'Order', 'Corynebacterium.like', NA, NA)
species <- except_treat(species, c('Clostridia bacterium', 'Lachnospiraceae oral'), 'Class', 'Phylum', 'Phylum', 'Clostridia', 'Firmicutes', 'Firmicutes')
species <- except_treat(species, 'Zetaproteobacteria bacterium', 'Class', 'Phylum', 'Phylum', 'Zetaproteobacteria', 'Proteobacteria', 'Proteobacteria')
species <- except_treat(species, 'Betaproteobacteria bacterium', 'Class', 'Phylum', 'Phylum', 'Betaproteobacteria', 'Proteobacteria', 'Proteobacteria')
species <- except_treat(species, 'Alphaproteobacteria bacterium', 'Class', 'Phylum', 'Phylum', 'Alphaproteobacteria', 'Proteobacteria', 'Proteobacteria')
species <- except_treat(species, 'Bacilli bacterium', 'Class', 'Phylum', 'Phylum', 'Bacilli', 'Firmicutes', 'Firmicutes')
# Traitement de 3 especes 'bacterium' ne pouvant etre fait avec la fonction except_treat()
ex <- which(species[, 'species']  %in% c('Acidobacteria bacterium', 'Actinobacteria bacterium', 'Tissierellia bacterium', 'Bacteroidetes bacterium', 'Tenericutes bacterium', 'Verrucomicrobia bacterium', 'Proteobacteria bacterium', 'Planctomycetes bacterium', 'Gammaproteobacteria bacterium', 'Firmicutes bacterium'))
species[ex, 'Phylum'] <- str_replace(species[ex, 'species'], '(.*) (.*)', '\\1')

all_species[, 1] <- species[, 1]

species %>%
  arrange(species) %>%
  identity -> species

all_species %>%
  arrange(species) %>%
  identity() -> all_species

#### Enregistrement de la dataframe slicee dans le fichier sliced_all_species_taxo.tsv et de 'species' dans taxo_species ####
write.table(all_species, "W:/Methode_Generic/output/Dataframe/sliced_all_species_taxo.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(species, "W:/Methode_Generic/output/Dataframe/taxo_species.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)