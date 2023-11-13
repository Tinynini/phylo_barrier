#library(tidyverse)
#library(tidytree)
#library(ape)

################################################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                                         #
# Utilite == generer les arbres taxonomiques a chaque niveaux par pruning de bac120_r95.tree   #
# Input == taxo_species.tsv, bac120_r95.tree et New_Parsed_taxonomy.tsv (V2)                   #
# Output == les 12 fichiers 'level_name[i]'.tree (V1) et 'level_name[i]'_version_alt.tree (V2) #
################################################################################################

#### Ouverture de bac120_r95.tree et de taxo_species.tsv & recuperation des donnees ####
tree <- read.tree('W:/Methode_Generic/data/bac120/bac120_r95.tree')
tree_df <- as_tibble(tree) # On passe au format tibble plus pratique a manipuler
# N.B. : Les labels de tips et de nodes se suivent sur une meme colonne lorsqu un arbre est au format tibble !!

species <- read_tsv('W:/Methode_Generic/output/Dataframe/taxo_species.tsv', col_types = 'cccccc') %>% 
  as.data.frame()

level_name <- unlist(colnames(species)) # On extrait les labels des 6 niveaux taxonomiques etudies pour pouvoir travailler a un niveau donne plus facilement

# On ajoute une colonne 'match' a species qui servira de marqueur lors du future join
match <- as.data.frame(matrix(data=1, nrow=nrow(species), ncol=1))
colnames(match) <- 'match'
species <- cbind(species, match) 

for (i in 1:6) # Permet de parcourir les 6 niveaux taxonomiques etudies (d espece a phylum)
{
  uni_level <- species[, c(i, 7)] # On extrait la colonne du niveau i et la colonne 'match'
  colnames(uni_level) <- c('level', 'match')
  
  #### Ouverture de New_Parsed_taxonomy.tsv & recuperation des donnees ####
  taxo <- read_tsv('W:/ninon-species/output/Table_taxonomie/New_Parsed_taxonomy.tsv', col_types = "cccccccc") %>%
    as.data.frame()
  
  small_taxo <- taxo[, c(7, i)] # On extrait les colonnes des labels d especes et du niveau i 
  
  #### Join de l arbre avec les colonnes de la table de taxonomie (== modification des labels de tips de l arbre) ####
  taxo_tree <- left_join(tree_df, small_taxo, by = c('label' = 'sseqid')) # On join l arbre et les colonnes de la table de taxonomie sur les colonnes de labels de part et d autre
  taxo_tree[1:Ntip(tree), 'label'] <- taxo_tree[1:Ntip(tree), colnames(taxo[i])] # On remplace le contenu de la colonne des labels par celui de la colonne du niveau i pour les labels de tips uniquement (surtout pas ceux de nodes !!)
  taxo_tree <- taxo_tree[, -5] # On supprime la colonne du niveau i
  Label <- taxo_tree[1:Ntip(tree),] # On extrait les lignes associers aux labels de tips
  
  Label %>% # Permet d isoler les lignes qui ont le meme contenu pour la colonne 'label' (== doublons) et de n en garder qu une seule a chaque fois
    arrange(label, parent) %>%
    group_by(label) %>%
    slice(1) -> labels # On cree une copie des les lignes isolees pour ne pas modifier l ensemble
  
  labels %>% # Permet de reordonner les lignes isolees de la meme facon que l ensemble
    arrange(node) %>%
    identity() -> labels
  
  # On extrait de nos lignes toutes celles qui ne font pas partie du groupe qu on a isole (== les doublons)
  tips <- which((Label[, 'node'] %>% pull()) %in% (labels[, 'node'] %>% pull()) == FALSE)
  new_tree <- as.phylo(taxo_tree) # On passe au format phylo pour pouvoir pruner l arbre
  new_tree <- drop.tip(new_tree, tips) # On prune l arbre en supprimant les lignes associees aux doublons
  taxo_tree <- as_tibble(new_tree) # On passe au format tibble plus pratique a manipuler
  
  #### Join de l arbre avec les colonnes du niveau i et de son partage (== transposition de l arbre a nos propres donnees) ####
  phylo_tree <- left_join(taxo_tree, uni_level, by = c('label' = 'level')) # On join l arbre et les colonnes du niveau i et de son partage sur la colonne des labels et celle du niveau i
  n_tips <- nrow(phylo_tree) - Nnode(new_tree) # Permet de recuperer le nombre de tips malgres l apparition de nombreux doublons
  phylo_tree <- unique(phylo_tree) # On dedoublonne notre arbre (operation impossible sans unifier d abbord les valeurs de partage)
  na_label <- is.na(phylo_tree[1:Ntip(new_tree), 'match']) # On refait le meme test en remettant a jour le nombre de tips
  level <- which(na_label == TRUE) # Cette fois on isole les lignes dont le partage n est pas renseigne (== celles qui ne matchent pas nos propres donnees)
  # N.B. : La presence de la colonne des partages du niveau i etait donc bien necessaire car elle seule nous permet de savoir quels tips ont ete matches ou non lors du join
  next_tree <- as.phylo(phylo_tree) # On passe au format phylo pour pouvoir pruner l arbre
  next_tree <- drop.tip(next_tree, level) # On prune l arbre en supprimant les lignes associees aux doublons
  phylo_tree <- as_tibble(next_tree) # On passe au format tibble plus pratique a manipuler
  
  #### Traitement special des labels de node pour les rendre uniques (donne lieu a une 2nd serie d arbres alternatifs)####
  phylo_tree %>% # On reordonne l arbre en fonction des labels parce que ca permet de separer automatiquement les labels de nodes deja uniques des autres et de regrouper ceux identiques
    arrange(label) %>%
    identity -> new_phylo_tree # On continue avec une copie de l arbre et non celui d origine pour avoir les 2 versions (avec ou sans traitement supplementaire pour les labels de nodes)
  
  nodes_exclus = length(grep("(.*)__(.*)", unlist(new_phylo_tree[, 'label']))) + 1 # On recupere le nombre de nodes deja uniques
  k <- 1 # Definit la valeur de la numerotation secondaire qu on va ajouter aux labels de node identiques
  
  for (j in nodes_exclus:Nnode(next_tree)) # On parcours les labels de nodes non uniques uniquement
  {
    if(new_phylo_tree[j + 1, 'label'] == new_phylo_tree[j, 'label']) # Si le label i est identique au label i + 1
    {
      new_phylo_tree[j, 'label'] <- str_glue("{new_phylo_tree[j, 'label']}_{k}") # On ajoute une numerotation secondaire k
      k <- k + 1 # Et on augmente k de 1 pour avancer dans la numerotation
    }
    else # Sinon
    {
      new_phylo_tree[j, 'label'] <- str_glue("{new_phylo_tree[j, 'label']}_{k}") # On ajoute quand même une numerotation secondaire k
      k <- 1 # Et on ramene k a la valeur 1 pour revenir au debut de la numerotation pour le groupe suivant de labels identiques
    }
  }
  
  new_phylo_tree %>% # On reordonne l arbre en fonction des numeros de node pour lui redonner sa structure d arbre
    arrange(node) %>%
    identity -> new_phylo_tree
  
  other_tree <- as.phylo(new_phylo_tree) # On repasse au format phylo pour pouvoir enregistrer l arbre sans risquer de l abimer
  
  #### Enregistrement des arbres ainsi obtenus dans des fichiers nominatifs ####
  path_start <- "W:/Methode_Generic/output/Arbre/"
  path_end <- ".tree"
  other_path_end <- "_version_alt.tree"
  # Les noms des fichiers sont definis par des variables
  file_name_1 <- str_glue("{path_start}{level_name[i]}{path_end}")
  file_name_2 <- str_glue("{path_start}{level_name[i]}{other_path_end}")
  
  write.tree(next_tree, file_name_1)
  write.tree(other_tree, file_name_2)
}