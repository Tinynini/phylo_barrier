#!/bin/bash

# INPUT
# species_ninon : listing des fichiers TSV associes aux especes traitees
# cluster : listing des genes retenus comme representants parmis ceux a traiter apres clusterisation

# OUTPUT
# centro : 1ere colonne de la matrice (= label de genes)
# species : 1ere ligne de la matrice ou de la fraction de matrice (= noms d'especes)
# share : corps de la matrice ou de la fraction de matrice

species_ninon='/home/ninon.robin/Methode_Generic/data/sort_species_ninon.txt'
cluster='/home/ninon.robin/Methode_Generic/data/clusters.txt'

centro='/home/ninon.robin/Methode_Generic/output/centro.tsv'
species='/home/ninon.robin/Methode_Generic/output/species.tsv'
share='/home/ninon.robin/Methode_Generic/output/share.tsv'

specie=$(cat ${species_ninon} | sed 's/^\/.*\/.*\/.*\//,/g' | sed 's/.tsv//g') 
echo -n -e 'SPECIES'${specie}\\\n | sed 's/ //g' > ${species}

centroid=$(cat ${cluster} | cut -f 10 | sort | uniq)

for centros in $(echo ${centroid})
do
    echo ${centros}
    echo ${centros} >> ${centro} # Ligne a commenter pour les fraction 2 a n si on travail sur des fractions !!
    query=$(cat ${cluster} | grep ${centros} | cut -f 9 | sort | uniq)
    query=$(echo ${query} | sed 's/ /|/g')

    for path in $(cat ${species_ninon})  
    do 
        match=$(cut -f 1 ${path} | grep -E -c ${query})
        if [ ${match} != 0 ]
        then 
            match=1
            echo "${path} has a match"
        fi

        echo -n ,${match} >> ${share}
    done 

    printf \\\n >> ${share}
done 