#library(tidyverse)

#########################################################################################
# Ninon ROBIN -- ninon.robin@inserm.fr                                                  #
# Utilite == generer les colonnes des partages inter-especes aux 6 niveaux taxonomiques #
# et les histogrammes des nombres d occurrence de valeurs de partages a chaque niveaux  #
# Input == taxo_species.tsv, uni_gene.tsv et les 6 matrices binaires                    #
# Output == Les 6 histogrammes (en FR et EN)                                            #
#########################################################################################

#### Ouverture de taxo_species.tsv et uni_gene.tsv & recuperation des donnees ####
species <- read_tsv('W:/Methode_Generic/output/Dataframe/taxo_species.tsv', col_types = 'cccccc') %>% 
  as.data.frame() 

uni_gene <- read_tsv('W:/Methode_Generic/output/Dataframe/uni_gene.tsv', col_types = 'c') %>% 
  as.data.frame()

#### Recuperation des partages au 6 niveaux taxonomiques etudies ####
level_name <- unlist(colnames(species)) # On extrait les labels des 6 niveaux taxonomiques etudies pour pouvoir travailler a un niveau donne plus facilement

m_path_start <- "W:/Methode_Generic/output/Matrice/Sliced_Matrix_"
m_path_end <- ".tsv"

level_shared <- as.data.frame(matrix(data=0, nrow=nrow(uni_gene), ncol=6))

for (i in 1:6)
{
  m_file_name <- str_glue("{m_path_start}{level_name[i]}{m_path_end}") # Le nom de fichier est definit par une variable
  matrix <- read.csv(file = m_file_name, header = TRUE, sep = "\t")
  
  level_share <- as.data.frame(matrix(data=0, nrow=nrow(uni_gene), ncol=1))
  
  for (j in 1:ncol(matrix))
  {
    level_share[j,] <- sum(matrix[, j])
  }
  
  level_shared[, i] <- level_share
}

# Netoyage de certains partages sortis a 0 au lieu de 1
for (i in 1:6)
{
  for (j in 1:nrow(level_shared))
  {
    if (level_shared[j, i] == 0)
    {
      level_shared[j, i] <- 1
    }
  }
}

#### Fonction pour generer les plots avec le titre et les labels en francais ####
generate_plot_fr <- function(shared_by, level_name)
{
  start <- 'Taxo_'
  end <- '_fr.png'
  title_start <- "Nombres d'occurrences des valeurs de partages inter-"
  path = "W:/Methode_Generic/output/Plot/Taxo_plot/FR"
  title <- str_glue("{title_start}{level_name}") # Le titre de l histogramme est definit par une variable
  ggplot(level_shared, aes(shared_by)) + geom_histogram(bins = max(shared_by)*2 - 1) + ggtitle(label = title) + xlab("Valeurs des partages") + ylab("Nombres d'occurences") 
  ggsave(str_glue("{start}{level_name}{end}"), plot = last_plot(), device = "png", path = path, width = 16, height = 8.47504)
}

#### Fonction pour generer les plots avec le titre et les labels en anglais ####
generate_plot_en <- function(shared_by, level_name)
{
  start <- 'Taxo_'
  end <- '_en.png'
  title_start <- "Inter-"
  title_end <- " sharing value occurences"
  path = "W:/Methode_Generic/output/Plot/Taxo_plot/EN"
  title <- str_glue("{title_start}{level_name}{title_end}") # Le titre de l histogramme est definit par une variable
  ggplot(level_shared, aes(shared_by)) + geom_histogram(bins = max(shared_by)*2 - 1) + ggtitle(label = title) + xlab("Sharing values") + ylab("Number of occurences") 
  ggsave(str_glue("{start}{level_name}{end}"), plot = last_plot(), device = "png", path = path, width = 16, height = 8.47504)
}

#### histogrammes des nombres d occurrences des valeurs de partage aux 6 niveaux taxonomiques etudies ####
for (i in 1:6)
{
  shared_by <- level_shared[, i]
  generate_plot_fr(shared_by, level_name[i]) # On lui applique la fonction generate_plot_fr()
  generate_plot_en(shared_by, level_name[i]) # On lui applique la fonction generate_plot_en()
}