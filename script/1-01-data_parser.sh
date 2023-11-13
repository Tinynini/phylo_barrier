#!/bin/bash

# INPUT
# select_genes : listing des sequences des genes a traiter
# res_path : chemin d'acces vers le dossier contenant les fichiers TSV obtenu en output de DIAMOND (un fichier par especes traitees)

select_genes='/home/ninon.robin/Methode_Generic/data/select_gene.ffn'

/opt/progs/vsearch --cluster_fast ${select_genes} --id 0.95 --uc /home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt

echo "VSEARCH CLUSTERING DONE"

cluster='/home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt'
cat ${cluster} | grep -E -v '^C' > /home/ninon.robin/Methode_Generic/data/clust.txt
  
IFS=$'\n'
    
for line in $(cat /home/ninon.robin/Methode_Generic/data/clust.txt)  
do 
    if [ -n "$(echo ${line} | grep -E '^S')" ]
    then 
        echo 'Processing centroid'

        start=$(echo ${line} | cut -f -9)
        echo -n -e ${start}\\\t >> /home/ninon.robin/Methode_Generic/data/clusters.txt
        echo ${line} | cut -f 9 >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    else
        echo 'Processing hit'

        echo ${line} >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    fi
done 

echo "Parsing clusters done"

res_path="/home/ninon.robin/diamond_ninon/"
find ${res_path} -name "*.tsv" > /home/ninon.robin/Methode_Generic/data/species_ninon.txt
species='/home/ninon.robin/Methode_Generic/data/species_ninon.txt'

i=0

while read path 
do 
    if [ $((i%100)) == 0 ];
    then 
        echo "${i}"; 
    fi

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