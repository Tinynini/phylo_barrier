#library(tidyverse)

###########################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                    #
# Utilite == generer les matrices d absence/presence GenexNiveau(xNiveau) #
# Input == Sliced_Matrix_species.tsv, taxo_species.tsv et uni_gene.tsv    #
# Output == 20 fichiers Sliced_Matrix_*.tsv                               #
###########################################################################

#### Ouverture de taxo_species.tsv et uni_gene.tsv & recuperation des donnees dans une dataframe ####
species <- read_tsv('W:/Methode_Generic/output/Dataframe/taxo_species.tsv', col_types = 'cccccc') %>% 
  as.data.frame()

level_name <- unlist(colnames(species)) # On extrait les labels des 6 niveaux taxonomiques etudies pour pouvoir travailler a un niveau donne plus facilement

uni_gene <- read_tsv('W:/Methode_Generic/output/Dataframe/uni_gene.tsv', col_types = 'c') %>% 
  as.data.frame()

uni_gene <- unlist(uni_gene)
n_gene <- length(uni_gene) 

#### Ouverture & pretraitrement de la matrice binaire GenexNiveau_i ####
for (i in 1:5) # Permet de parcourir les 5 niveaux taxonomiques etudies (d espece a classe)
{
  for (j in (i + 1):6) # Permet de parcourir en parrallele du niveau i les j niveaux suivants
  {
    curr_level <- as.data.frame(unique(species[, c(j, i)])) # On extrait simultanement les colonnes des niveaux j et i 
    
    if (i == 5) 
    { 
      na_level_2 <- which(is.na(curr_level[, 2]) == TRUE) # On isole les lignes du niveau i contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_2),] # On supprimr ces lignes
      
      #curr_level <- curr_level[-c(*ligne(s) a supprimer s il y en a*),]
      
      now_level <- as.data.frame(sort(unique(curr_level[, 2]))) # On extrait la colonne du niveau i 
      colnames(now_level) <- level_name[i]
    }
    
    else if (i == 4)
    {
      na_level_2 <- which(is.na(curr_level[, 2]) == TRUE) # On isole les lignes du niveau i contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_2),] # On supprimr ces lignes
  
      na_level_1 <- which(is.na(curr_level[, 1]) == TRUE) # On isole les lignes du niveau j contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_1),] # On supprime ces lignes
      
      curr_level <- curr_level[-c(31),] # /!\ les lignes a supprimer ici dependent des genes et especes traites !
      
      now_level <- as.data.frame(sort(unique(curr_level[, 2]))) # On extrait la colonne du niveau i 
      colnames(now_level) <- level_name[i]
    }
    
    else if (i == 3)
    {
      na_level_2 <- which(is.na(curr_level[, 2]) == TRUE) # On isole les lignes du niveau i contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_2),] # On supprimr ces lignes
      
      na_level_1 <- which(is.na(curr_level[, 1]) == TRUE) # On isole les lignes du niveau j contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_1),] # On supprime ces lignes
      
      curr_level <- curr_level[-c(72),] # /!\ les lignes a supprimer ici dependent des genes et especes traites !
      
      now_level <- as.data.frame(sort(unique(curr_level[, 2]))) # On extrait la colonne du niveau i 
      colnames(now_level) <- level_name[i]
    }
    
    else if (i == 2) 
    {
      na_level_2 <- which(is.na(curr_level[, 2]) == TRUE) # On isole les lignes du niveau i contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_2),] # On supprimr ces lignes
      
      now_level <- as.data.frame(sort(unique(curr_level[, 2]))) # On extrait la colonne du niveau i 
      colnames(now_level) <- level_name[i]
      
      na_level_1 <- which(is.na(curr_level[, 1]) == TRUE) # On isole les lignes du niveau j contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_1),] # On supprime ces lignes
    }
    
    else 
    {
      now_level <- as.data.frame(sort(unique(curr_level[, 2]))) # On extrait la colonne du niveau i 
      colnames(now_level) <- level_name[i]
      
      na_level_1 <- which(is.na(curr_level[, 1]) == TRUE) # On isole les lignes du niveau j contenant NA (== valeur non renseignee)
      curr_level <- curr_level[-c(na_level_1),] # On supprime ces lignes
    }
    
    curr_level %>% # On reordonne a present les 2 colonnes en fonction de celle du niveau j
      arrange(level_name[j]) %>%
      identity() -> curr_level
    
    uni_level <- unlist(as.data.frame(sort(unique(curr_level[, 1])))) # On extrait la colonne du niveau j
    n_level <- length(uni_level)

    # Obtention d une matrice formatee aux dimensions de curr_level
    path_start <- "W:/Methode_Generic/output/Matrice/Sliced_Matrix_"
    path_end <- ".tsv"
    file_name <- str_glue("{path_start}{level_name[i]}{path_end}") # Le nom de fichier est definit par une variable

    curr_matrix <- read.csv(file_name, header = TRUE, sep = "\t") # On ouvre la matrice depuis la liste de fichier

    temp_matrix <- cbind(now_level, curr_matrix) # On fusionne la colonne du niveau i a notre matrice
    
    rm(curr_matrix) # On supprime curr_matrix desormais inutile pour faire de la place

    temp_matrix <- left_join(curr_level, temp_matrix, by = NULL) # On join les colonnes des niveaux j et i et la matrice sur les colonnes du niveau i de part et d autre
    temp_matrix <- as.matrix(temp_matrix[, -c(1,2)]) # On supprime les colonnes des niveaux j et i

    #### Creaction d une matrice pseudo-binaire (0/1~n) d absence/presence des genes au niveau i au sein du niveau j ####
    cross_matrix <- matrix(data = 0, nrow = n_level, ncol = n_gene)

    colnames(cross_matrix) <- uni_gene
    rownames(cross_matrix) <- uni_level

    for (k in 1:n_level) # On parcourt les k representants distincts du niveau j
    {
      # On va les chercher dans la colonne triee et dedoublonnee qu on a extrait precedement
      to_set <- which(curr_level[, 1] %in% uni_level[k]) # On isole les occurrences du representant k au sein du bloc 'niveau j + niveau i'
      l <- length(to_set)
      mat <- temp_matrix[c(to_set),] # On extrait de notre matrice binaire les lignes de meme indice que les occurences trouvees
      # Pour la ligne associee au representant k dans la matrice pseudo-binaire :
      if (l > 1) # S il y a plus d une occurrence
      {
        cross_matrix[k,] <- colSums(mat) # On lui assigne comme contenu la somme des lignes extraites si-avant de la matrice binaire
      }

      else # Sinon
      {
        cross_matrix[k,] <- mat # On lui assigne comme contenu celui de l unique ligne extraite si-avant de la matrice binaire
      }
    } # N.B : Je sais c est un peu complique mais in fine ca donne un matrice Genex'Level j' avec 0 s il y a pas de match ou le nombre de representants du niveau i au sein du representant du niveau j se partageant le gene s il y a un match

    rm(temp_matrix) # On supprime temp_matrix desormais inutile pour faire de la place
   
    #### Enregistrement de la matrice pseudo-binaire ainsi obtenue dans un fichier nominatif ####
    next_file_name <- str_glue("{path_start}{level_name[i]}_{level_name[j]}{path_end}") # Le nom de fichier est definit par une variable
    write.table(cross_matrix, next_file_name, sep = '\t', row.names = TRUE, col.names = TRUE)

    #### Obtention & enregistrement de la matrice binaire au niveau j pour j = i + 1 dans un fichier nominatif ####
    if (j == i + 1)
    {
      for (o in 1:ncol(cross_matrix))
      {
        for (p in 1:nrow(cross_matrix))
        {
          if (cross_matrix[p, o] != 0)
          {
            cross_matrix[p, o] <- 1
          }
        }
      }

     new_file_name <- str_glue("{path_start}{level_name[j]}{path_end}") # Le nom de fichier est definit par une variable
     write.table(cross_matrix, new_file_name, sep = '\t', row.names = TRUE, col.names = TRUE)
    }

    rm(cross_matrix) # On supprime cross_matrix pour faire de la place a chaque fin de tour de la boucle j
  }
}