# Explication sur le contenu de ce github :

Ce github donne accès à l'ensemble des scripts R et bash permettant d'appliquer ma méthode d'analyse des partages de gènes bactériens de façon "casi" générique. La numérotation des scripts permet de savoir dans quel ordre ils doivent être utilisés, le scripts 2-00-main.R permettant de lancer directement l'ensemble des scripts R.

**Les scripts bash :**

Le 1er script a 3 tâches (+ une tâche optionelle) : 1) Clusteriser dans VSearch les séquences des gènes que l'on a choisi de traiter ; 2) Traiter le fichier output de ce clustering pour ne conserver que les parties utiles dans le bon formatage pour la suite ; 3) Supprimer les fichiers vides de la liste des fichiers TSV obtenus en output de DIAMOND (; 4) Fractionner la liste des fichiers TSV pour une obtention plus rapide de la matrice lorsque le nombre de gènes traités est élevé).

Le 2ème script sert à obtenir séparément le corps de la matrice gènesXespèces, ces noms de lignes (= gènes), et ces noms de colonnes (= espèces), soit pour l'ensemble des espèces, soit pour une fraction donnée de celles-ci. Ce même script peut donc être utiliser de 2 façons : 1) En le faisant tourner 1 fois sur la liste complète des fichiers TSV ; 2) En le faisant tourner n fois en parallèle sur les n fractions de cette liste. Le fichier centro.tsv dans lequel sont stoqués les noms de lignes (= gènes) ne doit être généré que pour la 1ère fraction. En pratique, cela signifit qu'il faut commenter la ligne de code associée pour les fractions de 2 à n. 

Le 3ème script sert à assembler la matrice à partir des outputs obtenus avec le script précédent. Il s'agit de la version à utiliser lorsque l'ensembles des espèces a été traiter d'un seul bloc.

Le 3ème script version bis sert à assembler la matrice à partir des outputs obtenus avec le script précédent. Il s'agit de la version à utilisr lorsque les espèces ont été traitées en n fractions

/!\ N.B. : 

La méthodes de fractionnement de la liste des fichiers TSV dans le 1er script et celle de la version bis du 3ème script doivent être adaptées en fonction du nombre n de fractions. Par ailleurs, dans le cas de la méthode de fractionnement, les numéros de lignes délimitant les fractions doivent s'enchainer. C'est à dire que le numéro de la ligne de départ d'une fraction i doit toujour être égale à 1 + celui de la ligne de fin de la fraction i - 1. Et le numéro de ligne de fin de la dernière fraction doit être égale au numéro de la dernière ligne non vide du fichier contenant la liste des fichiers TSV. 

**Les scripts R :**

L'ensemble des scripts R peut être lancer d'un seul coup via le script 2-00-main.R, mais ils peuvent aussi être lancer un par un, au grés de l'utilisateur.

2-01-02_alt.R : Sert à nétoyer la matrice et à récupérer les noms d'espèce dans une nouvelle structure.
	
2-03-Taxonomy_Parser.R : Sert à nétoyer les 2 versions de la table de taxonomie. 
--> /!\ Un traitement vise à éviter certains doublons lors du join dans le script suivant par l'utilisation successif d'une fonction sur différents cas d'espèces dont la taxonomie semble contenir une erreur ou bien être obsolète. De nouveaux cas similaires sont succeptibles de se présenter si des espèces que je n'ai encore jamais traitées devaient sortir en output de DIAMOND. Je vous invite à appliquer la fonction utilisée lors de ce traitement à ces nouveaux cas lorsque c'est possible pour limiter les cas à traiter déjà nombreux dans le script suivant.
	
2-04-Taxo_Join.R : Sert à associer leur taxonomies à un maximum de nos espèce dans la nouvelle structure (présence de traitements à adapter en fonction des données !).
	
2-05-Matrix_Preparator.R : Sert préparer l'obtention de l'ensemble des matrices binaires et pseudo-binaires dans le script suivant.
	
2-06-Matrix_Generator.R : Sert à obtenir les dites matrices (présence d'un traitement à adapter en fonction des données !)
	
2-07-Matrix_Plot.R : Sert à générer les barplots de présence/absence d'un gène donné à un niveau j au sein de chacun des ses représentant à un niveau inférieur i.
	
2-08-Taxo_Filter_Plot.R : Sert à générer les histogrammes des nombres d'occurrences des valeurs de partages de gènes.
	
2-09-Tree_Pruner.R : Sert à obtenir les arbres de nos espèces aux 6 niveaux taxonomique à partir de celui associé à la table de taxonomie.
		
2-10-Tree_Plot.R : Sert à obtenir aux 6 niveaux taxonomiques les sous-arbres associés à chaques partages et le plot de ceux-ci pour 4 gènes données, les histogrammes des distances totales des sous-arbres et les plots de l'arbre complet avec et sans les sous-arbres visibles. Les 4 gènes dont on plots les sous-abres doivent être choisis afin de représenter un partage important, 2 partages moyens distincts, et un petit partage.

/!\ N.B. :

Les scripts 2-07-Matrix_Plot.R et 2-10-Tree_Plot.R inclus tous 2 des plots de graphs de gènes données (et non sur tout les gènes). Pour choisir les gènes testés dans ces 2 cas, il faut se rendre directement dans les scripts pour modifier les gènes sélectionnés, d'autant que ceux choisis pour les données précédentes ne conviennent pas toujours pour les données suivantes (numéro d'emplacement dans la matrix précédente dépassant les dimmensions de la nouvelle par exemple, gènes trop ou trop peu partagés pour que les plots soient exploitables, etc.). Les instructions pour opérer ces modifications sans risque sont fournies dans le script 2-00-main.R et de façon plus précise directement dans les scripts concernés, aux niveaux-mêmes de ces plots.

Les scripts 2-04-Taxo_Join.R et 2-06-Matrix_Generator.R inclus tous 2 des traitements qui peuvent nécessités d'être adaptés car dépendant des gènes et espèces traités. Pour le script 2-06-Matrix_Generator.R, un seul traitement est concerné, et les modifications pouvant s'averer nécessaire touchent à l'exclusion de certaines lignes dans une structure temporaire pour certains niveaux taxonomiques. Pour l'autre script, il y a 4 traitements différents qui peuvent nécessités d'êtres adaptés, dont un directement au niveau d'une fonction. Les instructions pour opérer ces modifications sans risque sont fournis dans le script 2-00-main.R et de façon plus précise directement dans les scripts concernés, aux niveaux-mêmes de ces traitements.
