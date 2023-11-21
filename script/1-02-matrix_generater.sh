#!/bin/bash

# INPUT
# species_ninon : listing des fichiers TSV associes aux especes traitees
# cluster : output du clustering VSearch netoye (voir script 1-01)

# OUTPUT
# centro : 1ere colonne de la matrice (= label de genes)
# species : 1ere ligne de la matrice ou de la fraction de matrice (= noms d'especes)
# share : corps de la matrice ou de la fraction de matrice

species_ninon='/home/ninon.robin/Methode_Generic/data/sort_species_ninon.txt'
cluster='/home/ninon.robin/Methode_Generic/data/clusters.txt'

centro='/home/ninon.robin/Methode_Generic/output/centro.tsv'
species='/home/ninon.robin/Methode_Generic/output/species.tsv'
share='/home/ninon.robin/Methode_Generic/output/share.tsv'

# On recupere la liste des noms d'especes a partir de celle des chemins d'acces vers les fichiers TSV
specie=$(cat ${species_ninon} | sed 's/^\/.*\/.*\/.*\//,/g' | sed 's/.tsv//g') 
echo -n -e 'SPECIES'${specie}\\\n | sed 's/ //g' > ${species}

centroid=$(cat ${cluster} | cut -f 10 | sort | uniq) # On recupere la liste des centroids (10eme champs)

for centros in $(echo ${centroid}) # Pour chaque centroids
do
    # On enregistre le centroid dans une nouvelle ligne du fichier output prevu a cet effet
    echo ${centros} >> ${centro} # Ligne a commenter pour les fraction 2 a n si on travail sur des fractions !!
    
    query=$(cat ${cluster} | grep ${centros} | cut -f 9 | sort | uniq) # on preleve les genes associes au centroid
    query=$(echo ${query} | sed 's/ /|/g') # on formate la liste de genes en vue d'un grep sur cette liste

    for path in $(cat ${species_ninon}) # Pour chaque especes bacteriennes
    do 
        # On cherche les genes selectionnes ci-avant dans le 1er champs du fichier TSV associe 
        match=$(cut -f 1 ${path} | grep -E -c ${query}) # On ne recupere que le nombre de matchs
        
        if [ ${match} != 0 ] # Si le nombre de matchs est different de zero 
        then 
            match=1 # On lui donne la valeur 1 par defaut (binarisation de notre future matrice)
        fi

         # On recupere le nombre de matchs dans une nouvelle colonne du fichier output prevu a cet effet
        echo -n ,${match} >> ${share}
    done 

    printf \\\n >> ${share} # On passe a la ligne suivante (un centroid == une ligne)
done 
