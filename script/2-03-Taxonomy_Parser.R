#library(tidyverse)

##############################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                       #
# Utilite == pretraiter les 2 versions de la table de taxonomie en parallele #
# Input == bac120_taxonomy_r95.tsv (V1) et bac120_taxonomy_r95_new.tsv (V2)  #
# Output == Parsed_taxonomy.tsv (V1) et New_Parsed_taxonomy.tsv (V2)         #
##############################################################################

#### Ouverture de bac120_taxonomy_r95.tsv et de bac120_taxonomy_r95_new.tsv & recuperation des donnees ####
taxonomy_V1 <- read.csv('W:/Methode_Generic/data/bac120/bac120_taxonomy_r95.tsv', sep = ';', header = FALSE)
taxonomy_V2 <- read.csv('W:/Methode_Generic/data/bac120/bac120_taxonomy_r95_new.tsv', sep = ';', header = FALSE)

colnames(taxonomy_V1) <- c('Occurrence', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
colnames(taxonomy_V2) <- c('sseqid', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

#### Modification de la nomenclature du join avec nos propre donnees ####
nomenclature_cleaner <- function(taxonomy)
{
  for (i in 2:ncol(taxonomy)) # Suppression des prefix indiquant le niveau taxonomique (superflus car deja indique par les noms de colonnes)
  {
    taxonomy[,i] <- str_replace(taxonomy[,i], pattern = '(...)(.*)', replacement = "\\2")
  }
  
  for (j in 1:nrow(taxonomy)) # Suppression des supplements de matricule en suffix chez les especes 'sp' (egalement superflus)
  {
    if (str_replace(taxonomy[j,'Species'], pattern = "(.*) (..)(.*)", replacement = "\\2") == 'sp')
    {
      taxonomy[j,'Species'] <- str_replace(taxonomy[j,'Species'], pattern = "(.*) (..)(.*)", replacement = "\\1\\ \\2")
    }
  }
  
  taxonomy %>% 
    arrange(Species) %>% 
    identity() -> taxonomy
}

#### Elimination de certaines sous-classifications propres a l'etude dont est tiree cette table (et parasites pour la notre !) ####
sub_class_cleaner <- function(taxonomy)
{
  taxonomy <- rev(unique(taxonomy))
  
  for (k in 2:6) # On commence a 2 car la colonne des especes n a pas besoin d etre traitee 
  {
    sub <- grep('(.*)_(.*)', taxonomy[, k])
    taxonomy[sub, k] <- str_replace(taxonomy[sub, k], '(.*)_(.*)', '\\1')
  }
  
  sub_espece1 <- grep('(.*)_(.*) (.*)', taxonomy[, 1]) # Cette fois on ne traite que la colonne des especes
  taxonomy[sub_espece1, 1] <- str_replace(taxonomy[sub_espece1, 1], '(.*)_(.*) (.*)', '\\1\\ \\3')
  sub_espece2 <- grep('(.*) (.*)_(.*)', taxonomy[, 1]) # Cette fois on ne traite que la colonne des especes
  taxonomy[sub_espece2, 1] <- str_replace(taxonomy[sub_espece2, 1], '(.*) (.*)_(.*)', '\\1\\ \\2')
  
  taxonomy <- unique(taxonomy)
}

#### Traitement preventif pour eviter la creation de certains doublons lors du join ####
sub_grep_level <- function(taxonomy, pattern, level, replacement_1, replacement_2)
{
  sub <- grep(pattern, taxonomy[, level])
  taxonomy[sub, 'Family'] <- replacement_1
  taxonomy[sub, 'Order'] <- replacement_2
  
  return(taxonomy)
}

prev_doublon_cleaner <- function(taxonomy) # Permet d'appliquer plusieurs fois sub_greb_level() avec differents parametrage
{
  taxonomy <- sub_grep_level(taxonomy, 'Bacillus', 'Genus', 'Bacillaceae', 'Bacillales')
  taxonomy <- sub_grep_level(taxonomy, 'Ruminococcus sp', 'Species', 'Ruminococcaceae', 'Oscillospirales')
  taxonomy <- sub_grep_level(taxonomy, 'Clostridium sp', 'Species', 'Clostridiaceae', 'Clostridiales')
  taxonomy <- sub_grep_level(taxonomy, 'Eubacterium sp', 'Species', 'Lachnospiraceae', 'Lachnospirales')
  taxonomy <- sub_grep_level(taxonomy, 'Leptolyngbya ohadii', 'Species', 'Leptolyngbyaceae', 'Leptolyngbyales')
  
  taxonomy <- unique(taxonomy)
}

#### Main ####
taxonomy_V1 <- nomenclature_cleaner(taxonomy_V1)
taxonomy_V2 <- nomenclature_cleaner(taxonomy_V2)
taxonomy_V1 <- taxonomy_V1[, -c(1, 2)]
taxonomy_V2 <- taxonomy_V2[, -c(2)]
taxonomy_V1 <- sub_class_cleaner(taxonomy_V1)
taxonomy_V2 <- sub_class_cleaner(taxonomy_V2)
taxonomy_V1 <- prev_doublon_cleaner(taxonomy_V1)
taxonomy_V2 <- prev_doublon_cleaner(taxonomy_V2)

#### Enregistrement des 2 versions de la table de taxonomie dans les fichier Parsed_taxonomy.tsv et New_Parsed_taxonomy.tsv ####
write.table(taxonomy_V1, "W:/Methode_Generic/output/Table_taxonomie/Parsed_taxonomy.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(taxonomy_V2, "W:/Methode_Generic/output/Table_taxonomie/New_Parsed_taxonomy.tsv", sep = '\t', row.names = FALSE, col.names = TRUE)