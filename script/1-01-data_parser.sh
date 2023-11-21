#!/bin/bash

# INPUT
# select_genes : listing des sequences des genes a traiter
# res_path : chemin d'acces vers le dossier contenant les fichiers TSV obtenu en output de DIAMOND (un fichier par especes traitees)

#### Clusterisation des genes traites avec VSearch ####
# /!\ Il faut utiliser le CDC mothur pour y avoir acces !!

select_genes='/home/ninon.robin/Methode_Generic/data/select_gene.ffn'

/opt/progs/vsearch --cluster_fast ${select_genes} --id 0.95 --uc /home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt

echo "VSEARCH CLUSTERING DONE"

#### Traitement de fichier obtenu en output du clustering ####

cluster='/home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt'
# On supprime les lignes qui listent les clusters en fin de fichier (celles qui commence par "C")
cat ${cluster} | grep -E -v '^C' > /home/ninon.robin/Methode_Generic/data/clust.txt
  
IFS=$'\n'
    
for line in $(cat /home/ninon.robin/Methode_Generic/data/clust.txt)  
do 
    if [ -n "$(echo ${line} | grep -E '^S')" ] # Pour les lignes associees a un centroid (qui commence par "S")
    then 
        echo 'Processing centroid'
        # On supprime le 10eme champs et on recopie le 9eme champs a la place
        start=$(echo ${line} | cut -f -9) 
        echo -n -e ${start}\\\t >> /home/ninon.robin/Methode_Generic/data/clusters.txt
        echo ${line} | cut -f 9 >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    else # Pour les ligne associees a un gene lambda
        echo 'Processing hit'
        # On ne fait aucun traitement particulie
        echo ${line} >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    fi
done 

echo "Parsing clusters done"

#### Suppression des fichier TSV vides ####

res_path="/home/ninon.robin/diamond_ninon/"

# On recupere les chemins d'acces de tous les fichiers TSV que contient res_path
find ${res_path} -name "*.tsv" > /home/ninon.robin/Methode_Generic/data/species_ninon.txt
species='/home/ninon.robin/Methode_Generic/data/species_ninon.txt'

i=0

while read path 
do 
    if [ $((i%100)) == 0 ]; # On ne print cet indicateur que toutes les 100 lignes (sinon ça fait trop)
    then 
        echo "${i}"; 
    fi
    # Si le fichier associe a la ligne traitee ne commence par par une ligne vide on le garde
    if [ -n "$(head -n 1 ${path})" ]
    then 
        echo ${path} >> /home/ninon.robin/Methode_Generic/data/sort_species_ninon.txt
    fi

    i=$(($i+1))
done < "${species}"

#### Methode pour fractionner la liste de fichier TSV si le nombre de genes traites est eleves ####

#du -a $(cat /home/ninon.robin/Data_recup_reduc/data/sort_species_ninon.txt) | sort -r -g > /home/ninon.robin/Data_recup_reduc/data/sort_DU.txt
#cat Data_recup_reduc/data/sort_DU.txt | sed -n '1,*num_ligne_fin_de_fraction*p' | cut -f 2 > Data_recup_reduc/data/Parsed_species/*nom_fraction_1*.txt
#... (<- Dupliquer ici la ligne precedente n-2 fois)
#cat Data_recup_reduc/data/sort_DU.txt | sed -n '*num_ligne_debut_de_fraction*,*num_derniere_ligne*p' | cut -f 2 > Data_recup_reduc/data/Parsed_species/*nom_fraction_n*.txt
