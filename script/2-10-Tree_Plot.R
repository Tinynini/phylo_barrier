#library(tidyverse)
#library(tidytree)
#library(ape)
#library(ggtree)

##############################################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                                       #
# Utilite == generer les sous-arbres d un gene donne, les histogrammes des distances totales #
# des sous-rabres et les arbres complets avec sous-arbres colorises dessus à chaque niveaux  #
# Input == taxo_species.tsv, les 12 fichiers 'level_name[i]'*.tree                           #
# et les 6 fichiers Sliced_matrix_'level_name[i]'.tsv                                        #
# Output == les 6 histogrammes (en FR et EN), les 6 arbres en version normale (en FR et EN)  #
# et en version colorisee (en FR et EN) et les sous-arbres (6 ou moins selon le gene teste)  #
##############################################################################################

#### Ouverture de taxo_species.tsv & recuperation des donnees ####
species <- read_tsv('W:/Methode_Generic/output/Dataframe/taxo_species.tsv', col_types = 'cccccc') %>% 
  as.data.frame() 

# On est oblige de modifier la nomenclature des noms d especes parce qu une modification automatique se fait au niveau des labels de tips de l arbre des especes 
species[, 'species'] <- str_replace(species[, 'species'], '(.*) (.*)', '\\1\\_\\2')
level_name <- unlist(colnames(species)) # On extrait les labels des 6 niveaux taxonomiques etudies pour pouvoir travailler a un niveau donne plus facilement

#### Preparation des futures listes dans lesquels seront reunies celles obtenues aux 6 niveaux taxonomiques ####
liste_uni_gene <- vector(mode = 'list', length = 6) # On prepare une liste des listes des genes et des distances totales de leurs sous-arbres aux 6 niveaux taxonomiques
liste_tree_liste <- vector(mode = 'list', length = 6) # On prepare une liste des listes des sous-arbres aux 6 niveaux taxonomiques
min_length <- vector(mode = 'list', length = 6) # On prepare une liste des distances totales de sous-arbres minimales aux 6 niveaux taxonomiques
max_length <- vector(mode = 'list', length = 6) # On prepare une liste des distances totales de sous-arbres maximales aux 6 niveaux taxonomiques

#### Fonction servant a la creation de nouvelles listes des sous-arbres par genes et de leurs distances totales ####
liste_generator <- function(tree, tibble_tree) # Il faut l arbre sous forme phylo et sous forme tibble en entree 
{
  n_Gene <- ncol(tibble_tree)
  trees <- vector(mode = 'list', length = n_gene) # On prepare une liste des sous-arbres 
  length <- as.data.frame(matrix(data = 0, nrow = n_gene, ncol = 1)) # On prepare une colonne des distances des sous-arbres  
  uni_gene <- cbind(uni_gene, length) # On ajoute cette colonne a celle des genes
  colnames(uni_gene) <- c('gene', 'length')
  
  l <- 1
  
  for (k in 5:n_Gene) # Permet de parcourir les k colonnes associees aux genes dans tibbled_tree (celles issues de la matrice)
  { # N.B. : On est donc oblige de demarrer a partir de la 5eme colonnes (les 4 1ere etant celles propres a l arbre)
    wanted_gene <- colnames(tibble_tree[, k]) # On recuppere le nom de le gene associe a la colonne k
    wanted_tip <- tibble_tree$label[tibble_tree[wanted_gene] == 1] # On recupere les labels de tips se partagent le gene (== les lignes pour lesquelles il y a "1" dans la colonne du gene)
    wanted_tip <- na.omit(wanted_tip) # On doit exclures les 'NA' qui sont apparement consideres par defauts comme correspondant au '1' recherche ci-dessus (ils correspondent aux lignes des labels de nodes dont on en veut surtout pas !)
    tree_gene <- keep.tip(tree, tip = wanted_tip) # On prune l arbre complet pour ne garder que les tips selectionnes ci-avant
    length <- sum(tree_gene$edge.length) # On somme les distances des branches du sous-arbre pour recuperer sa distance totale
    uni_gene[k - 4, 'length'] <- length # On la stock dans la nouvelle colonne de uni_gene
    # N.B. : Les genes sont ordonnes de la meme facon dans tibble_tree et uni_gene donc il suffit de parcourir uni_gene en parallele (en partant bien de 1 et non plus de 5 cette fois) pour etre toujour a la bonne ligne
    trees[[l]] <- tree_gene # On stock le sous-arbre dans la liste des sous-arbre
    l <- l + 1
  } # Comme return ne peut s appliquer qu a une seule variable on est oblige de stocker temporairement les 2 listes ensemble
  
  liste <- list(trees, uni_gene) 
  return(liste)
}

#### Fonction servant a suppimer les sous-arbres vides dans une liste de sous_arbres ####
liste_parser <- function(trees, uni_gene) # Il faut la liste des sous-arbre et uni_gene (le block : gene + distance) en entree
{
  n_gene <- nrow(uni_gene)
  tree_list <- vector(mode = 'list', length = n_gene) # On prepare une nouvelle liste des sous-arbres 
  
  l <- 1
  
  for (k in 1:length(trees)) # Permet de parcourir les k sous-arbres de la liste
  {
    if (is.null(trees[[k]]) == FALSE) # Si le sous-arbre k n est pas vide
    {
      tree_list[[l]] <- trees[[k]] # On le copie dans la nouvelle liste
      l <- l + 1
    }
  } # On renomme les sous_arbres en fonction des genes auxquels ils sont associes (sinon il faudrait se referer constamment a uni_gene pour savoir a quel gene est associe un sous-arbre)
  
  names(tree_list) <- uni_gene[, 'gene'] 
  return(tree_list)
}

#### Main ####
for (i in 1:6) # Permet de parcourir les 6 niveaux taxonomiques etudies (d espece a phylum)
{
  #### Ouverture des arbres du niveau i depuis leurs fichiers nominatifs & preparation des donnees ####
  path_start <- "W:/Methode_Generic/output/Arbre/"
  path_end <- ".tree"
  other_path_end = "_version_alt.tree"
  # Les noms des fichiers sont definis par des variables
  file_name_1 <- str_glue("{path_start}{level_name[i]}{path_end}") 
  file_name_2 <- str_glue("{path_start}{level_name[i]}{other_path_end}")
  tree <- read.tree(file_name_1) # Arbre sans le traitement supplementaire des labels de nodes
  other_tree <- read.tree(file_name_2) # Arbre avec le traitement supplementaire des labels de nodes
  
  ggtree(tree) + ggtitle(str_glue("Arbre des {level_name[i]}"))
  ggsave(str_glue("tree_{level_name[i]}_fr.png"), plot = last_plot(), device = "png", path = "W:/Methode_Generic/output/Plot/Tree_plot/Arbres/FR", width = 16, height = 8.47504)
  ggtree(tree) + ggtitle(str_glue("{level_name[i]} tree"))
  ggsave(str_glue("tree_{level_name[i]}_en.png"), plot = last_plot(), device = "png", path = "W:/Methode_Generic/output/Plot/Tree_plot/Arbres/EN", width = 16, height = 8.47504)
  
  tibble_tree <- as_tibble(tree) # On passe au format tibble plus pratique a manipuler
  other_tibble_tree <- as_tibble(other_tree) # On passe au format tibble plus pratique a manipuler

  ### Ouverture & traitement de la matrice binaire associee au niveau i et de uni_level.tsv ####
  m_path_start <- "W:/Methode_Generic/output/Matrice/Sliced_Matrix_"
  m_path_end <- ".tsv"
  m_file_name <- str_glue("{m_path_start}{level_name[i]}{m_path_end}") # Le nom de fichier est definit par une variable
  gene_matrix <- read.csv(file = m_file_name, header = TRUE, sep = "\t") # On ouvre la matrice depuis la liste de fichier

  uni_gene <- read_tsv('W:/Methode_Generic/output/Dataframe/uni_gene.tsv', col_types = 'c') %>%
    as.data.frame()

  n_gene <- nrow(uni_gene)

  uni_level <- as.data.frame(rownames(gene_matrix)) # On extrait la colonne du niveau i

  colnames(uni_level) <- level_name[i]
  n_level <- nrow(uni_level)

  #### Join de l arbre et de la matrice & preparation de nouvelles listes ####
  gene_matrix <- cbind(uni_level, gene_matrix) # On combine la colonne du niveau i a la matrice en vue du join avec l arbre du niveau i

  if (i == 1)
  {
    gene_matrix[, 'species'] <- str_replace(gene_matrix[, 'species'], '(.*) (.*)', '\\1\\_\\2')
  }

  tibble_tree <- left_join(tibble_tree, gene_matrix, by = c('label' = level_name[i])) # On join la matrice au 1er arbre sur les colonnes du niveau i et des labels
  other_tibble_tree <- left_join(other_tibble_tree, gene_matrix, by = c('label' = level_name[i])) # On join la matrice au 2nd arbre sur les colonnes du niveau i et des labels

  #### Creation des listes des sous-arbres et de leurs distances totales par genes ####
  liste <- liste_generator(tree, tibble_tree) # On genere la listes de sous_arbres et la nouvelle colonne d uni_gene de leurs distances totales pour le 1er arbre
  other_liste <- liste_generator(other_tree, other_tibble_tree) # Idem pour le 2nd arbre
  # On recupere separement la liste des sous arbre et uni_gene pour les 2 arbres
  trees <- liste[[1]]
  other_trees <- other_liste[[1]]
  uni_gene <- liste[[2]]
  other_uni_gene <- other_liste[[2]]

  #### Exemple de plot des sous_arbres de genes tres, moyennement (x2) et peu partages (pour le 1er arbre uniquement parce que c est pareil si on le fait avec l autre) ####
  # Pour definir les noms et destinations de fichiers pour l enregistrement
  debu_b <- "W:/Methode_Generic/output/Plot/Tree_plot/Sous_arbres/rep~big_share/Sub_tree_"
  debu_m_1 <- "W:/Methode_Generic/output/Plot/Tree_plot/Sous_arbres/rep~medium_share-1/Sub_tree_"
  debu_m_2 <- "W:/Methode_Generic/output/Plot/Tree_plot/Sous_arbres/rep~medium_share-2/Sub_tree_"
  debu_s <- "W:/Methode_Generic/output/Plot/Tree_plot/Sous_arbres/rep~small_share/Sub_tree_"
  fine <- ".png"
  
  #### Fonction servant a effectuer le plot des sous-arbres d un gene donne ####
  sub_tree_plot <- function(gene, debu) # Il faut la liste des sous-arbre et uni_gene (le block : gene + distance) en entree
  { # N.B. : C est normal* qu a partir du niveau des genres les sous-arbres puissent etre vides pour les petits partages
    png(str_glue("{debu}{level_name[i]}{fine}"), height = 1017, width = 1920, pointsize = 20)
    plot.phylo(trees[[gene]], show.node.label = TRUE, main = uni_gene[gene, 1], sub = uni_gene[gene, 2])
    dev.off()
  } # */!\ Sous-arbre vide <=> pas d embranchement <=> pas de partage <=> partage = 1 (si patage > 1 alors gros problemo !!)
  
  # N.B. : Pour chaque categorie il suffit de modifier la 1ere variable pour tester un autre gene comme representant
  sub_tree_plot(21, debu_b) # Gros partage (grosse distance totale dans uni_gene au niveau des phyla)
  sub_tree_plot(10, debu_m_1) # Partage moyen 1 (distance totale mediane dans uni_gene au niveau des phyla)
  sub_tree_plot(6, debu_m_2) # Partage moyen 2 (distance totale mediane dans uni_gene au niveau des phyla)
  sub_tree_plot(4, debu_s) # Petit partage (petite distance totale dans uni_gene au niveau des phyla)
  
  #### Suppresion des sous_arbres vides et de leurs distances totales (genant pour la suite) ####
  err <- which(uni_gene[, 'length'] == 0.000) # On isole les lignes associees a des distances totales null (celles des sous-arbres vides) pour le 1er arbre
  other_err <- which(other_uni_gene[, 'length'] == 0.000) # Idem pour le 2nd arbre
  uni_gene <- uni_gene[-c(err),] # On supprime ces lignes de uni_gene pour le 1er arbre
  other_uni_gene <- other_uni_gene[-c(other_err),] # Idem pour le 2nd arbre
  tree_list <- liste_parser(trees, uni_gene) # On genere la nouvelle liste des sous-arbres sans ceux vides pour le 1er arbre
  other_tree_list <- liste_parser(other_trees, other_uni_gene) # Idem pour le 2nd arbre

  #### Histogrammes des distances totales des sous-arbres (la encore c est identique pour les 2 arbres donc on le fait que pour le 1er) ####
  level_length <- uni_gene['length']
  # Pour definir les noms et destinations de fichiers pour l enregistrement
  start <- 'Dist_'
  end_fr <- '_fr.png'
  end_en <- '_en.png'
  path_fr = "W:/Methode_Generic/output/Plot/Distance_plot/FR"
  path_en = "W:/Methode_Generic/output/Plot/Distance_plot/EN"
  # Pour definir les titres de plots
  title_fr <- "Nombres d'occurrences des valeurs de distances inter-"
  title_start_en <- "Inter-"
  title_end_en <- " distance value occurences"
  # On fait un premier plot avec le titre et les legendes en francais puis un second avec le titre et les legendes en anglais
  ggplot(level_length, aes(length)) + geom_histogram(bins = n_gene) + ggtitle(label = str_glue("{title_fr}{level_name[i]}")) + xlab("Valeurs des distances") + ylab("Nombres d'occurrences")
  ggsave(str_glue("{start}{level_name[i]}{end_fr}"), plot = last_plot(), device = "png", path = path_fr, width = 16, height = 8.47504)
  ggplot(level_length, aes(length)) + geom_histogram(bins = n_gene) + ggtitle(label = str_glue("{title_start_en}{level_name[i]}{title_end_en}")) + xlab("Distances values") + ylab("Number of occurences")
  ggsave(str_glue("{start}{level_name[i]}{end_en}"), plot = last_plot(), device = "png", path = path_en, width = 16, height = 8.47504)

  #### Plot des sous-arbres des especes par genes sur l arbre complet (Pour le 2nd arbre cette fois parce que ca ne peut pas fonctionner sans le traitement supplementaire des labels de nodes !!) ####
  liste <- vector(mode = 'list', length = length(other_tree_list)) # On prepare une nouvelle liste
  for (m in 1:length(other_tree_list)) # Permet de parcourir les m sous-arbre de la liste
  {
    wanted_tree <- as_tibble(other_tree_list[[m]]) # On recupere le sous-arbre m sous la forme d un tibble
    root <- which(is.na(wanted_tree['branch.length']) == TRUE) # On isole la ligne associee a sa racine dont la distance est la seule non renseignee (== 'NA')
    label <- which(other_tibble_tree$label %in% wanted_tree[root, 'label']) # On isole les lignes de meme label que sa racine dans l arbre complet
    liste[m] <- other_tibble_tree[label, 'node'] # On recuppere les numeros de node associes a ces ligne dans la nouvelle liste
  }
  liste <- t(as.data.frame(unique(liste))) # On transforme la liste ainsi remplie en dataframe dedoublonnee (ca necessite une transposition)
  type <- as.data.frame(matrix(data = 1:length(liste), nrow = length(liste), ncol = 1)) # On prepare une nouvelle colonne remplie avec des nombres allant de 1 au nombre de sous-arbres
  # N.B. : La colonne 'type' sert a donner des types distinct aux sous-arbres lors du plot via une numerotation pour pouvoir les coloriser tous differement sur l arbre complet
  liste <- cbind(liste, type) # On ajoute cette colonne a notre dataframe dedoublonnee
  names(liste) <- c('node', 'type')

  # Pour definir les noms et destinations de fichiers pour l enregistrement
  debut <- 'Tree_'
  fin_fr <- '_fr.png'
  fin_en <- '_en.png'
  dir_fr = "W:/Methode_Generic/output/Plot/Tree_plot/Arbre_sous_arbres/FR"
  dir_en = "W:/Methode_Generic/output/Plot/Tree_plot/Arbre_sous_arbres/EN"
  # Pour definir les titres de plots
  titre_deb_fr <- "Sous-arbres "
  titre_fin_fr <- "/ARG"
  titre_fin_en <- "/ARG sub-trees"
  # On fait un premier plot avec le titre en francais puis un second avec le titre en anglais
  # N.B. : geom_highlight permet de coloriser les sous-arbres en fonction du type associe. Il fallait donc definir autant de types differents qu il y a de sous-arbres pour attribuer une teinte unique a chacun
  ggtree(other_tree) + geom_hilight(data = liste, mapping = aes(node = node, fill = type)) + ggtitle(str_glue("{titre_deb_fr}{level_name[i]}{titre_fin_fr}"))
  ggsave(str_glue("{debut}{level_name[i]}{fin_fr}"), plot = last_plot(), device = "png", path = dir_fr, width = 16, height = 8.47504)
  ggtree(other_tree) + geom_hilight(data = liste, mapping = aes(node = node, fill = type)) + ggtitle(str_glue("{level_name[i]}{titre_fin_en}"))
  ggsave(str_glue("{debut}{level_name[i]}{fin_en}"), plot = last_plot(), device = "png", path = dir_en, width = 16, height = 8.47504)

  #### Stockages des listes obtenues au niveau taxonomique i dans les listes prevues a cet effet ####
  liste_uni_gene[[i]] <- uni_gene # On stock uni_gene dans la liste prevue pour ca
  liste_tree_liste[[i]] <- tree_list # On stock la liste des sous-arbre dans la liste prevue pour ca
  min_length[[i]] <- min(uni_gene[, 2]) # On recupere la valeur de distance totale minimale dans la liste prevue pour ca
  max_length[[i]] <- max(uni_gene[, 2]) # On recupere la valeur de distance totale maximale dans la liste prevue pour ca
}

# On enregistre les listes de sous-arbres et de genes dans un fichier RData
save(liste_uni_gene,  liste_tree_liste, file = "W:/Methode_Generic/output/listes.RData")