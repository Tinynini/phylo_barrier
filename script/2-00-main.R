library(tidyverse)
library(tidytree)
library(ape)
library(ggtree) 

#### Main : ####

source('W:/Methode_Generic/script/2-01-02_alt.R') 
source('W:/Methode_Generic/script/2-03-Taxonomy_parser.R')

source('W:/Methode_Generic/script/2-04-Taxo_Join.R') # /!\ 4 traitements peuvent necessiter d etre adaptes 
# en fonction des especes traites et pour 2 d entre eux de la taxonomie utilisee en ref !

source('W:/Methode_Generic/script/2-05-Matrix_Preparator.R')

source('W:/Methode_Generic/script/2-06_Matrix_Generator.R') # /!\ Les 3 1ers "if" (i = 3, 4 et 5) peuvent 
# necessiter d etre adaptes en fonction des genes et des especes traites !

source('W:/Methode_Generic/script/2-07-Matrix_Plot.R') # Il faut aller dans le script pour modifier le gene
# selectionne (celui dont on genere les barplots)

source('W:/Methode_Generic/script/2-08-Taxo_Filter_Plot.R') 
source('W:/Methode_Generic/script/2-09-Tree_Pruner.R')

source('W:/Methode_Generic/script/2-10-Tree_Plot.R') # Pour le plot des sous-arbres, il faut aller dans le 
# script pour modifier les representants slectionnes pour chaques categories de partage.

# N.B. : Les traitements pouvant necessiter une modification dans les scripts 04 et 06 sont signalises
# et leur fonctionnement detaillees dans les commentaires. Les plots s appliquant a un seul gene dans les
# scripts 07 et 10 sont signalises et les procedures pour modifier les genes selectionnes bien specifiees.
# Je n ai pas pu effectuer les modifications qui permettrait de selectionner ces genes lorsqu on lance les
# scripts, mais j ai simplifie au maximum ces procedures pour compenser la gene occasionnee.

# assign(x, value, pos = -1, envir = as.environment(pos), inherits = FALSE, immediate = TRUE) Solluce ??