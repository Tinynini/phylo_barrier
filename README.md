# Explication sur le contenu de ce github :

Ce github donne accès à l'ensemble des scripts R et bash permettant d'appliquer ma méthode d'analyse des partages de gènes bacterien de façon "casi" générique. La numérotation des scripts permet de savoir dans quel ordre ils doivent être utilisés, le scripts 2-00-main.R permettant de lancer directement l'ensemble des scripts R.

**Les scripts bash :**

Le 1er script a 3 tâches (+ une tâche optionelle) : 1) Clusteriser dans VSearch les séquences des gènes que l'on a choisi de traiter ; 2) Traiter le fichier output de ce clustering pour ne conserver que les parties utiles dans le bon formatage pour la suite ; 3) Supprimer les fichiers vides de la liste des fichiers de clustering obtenus en output de DIAMOND (; 4) Fractionner la liste des fichiers de clustering pour une obtention plus rapide de la matrice lorsque le nombre de gènes traités est élevé).

Le 2ème script sert à obtenir le corps de la matrice gènesXespèces avec ces noms de lignes, soit pour l'ensemble des espèces, soit pour une fraction donnée de celles-ci. Ce même script peut donc être utiliser de 2 façons : 1) En le faisant tourner 1 fois sur la liste complètes des fichiers de clustering ; 2) En le faisant tourner n fois en parallèle sur les n fractions de cette liste. Le fichier centro.tsv dans lequel est stoquée la 1ère colonne de la matrice, celle des noms de lignes (= gènes), ne doit être généré que pour la 1ère fraction. En pratique cela signifit qu'il faut commenter la ligne de code ... A FINIR

A FAIRE

Le 3ème script

Le 3ème script (version bis)

**Les scripts R :**
