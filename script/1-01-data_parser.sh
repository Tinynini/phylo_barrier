#!/bin/bash

/opt/progs/vsearch --cluster_fast /home/ninon.robin/Methode_Generic/data/select_gene.ffn --id 0.95 --uc /home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt

cluster='/home/ninon.robin/Methode_Generic/data/cluster_fast_select_gene_0.95.txt'
cat ${cluster} | grep -E -v '^C' > /home/ninon.robin/Methode_Generic/data/clust.txt
  
IFS=$'\n'
    
for line in $(cat /home/ninon.robin/Methode_Generic/data/clust.txt)  
do 
    if [ -n "$(echo ${line} | grep -E '^S')" ]
    then 
        start=$(echo ${line} | cut -f -9)
        echo -n -e ${start}\\\t >> /home/ninon.robin/Methode_Generic/data/clusters.txt
        echo ${line} | cut -f 9 >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    else
        echo ${line} >> /home/ninon.robin/Methode_Generic/data/clusters.txt
    fi
done 

res_ninon_path="/home/ninon.robin/diamond_ninon/"
find ${res_ninon_path} -name "*.tsv" > /home/ninon.robin/Methode_Generic/data/species_ninon.txt
cat /home/ninon.robin/Methode_Generic/data/species_ninon.txt | sort -f -t '/' -k 4 > /home/ninon.robin/Methode_Generic/data/full_species_ninon.txt

full_species='/home/ninon.robin/Methode_Generic/data/full_species_ninon.txt'

while read path 
do 
    if [ -n "$(head -n 1 ${path})" ]
    then 
        echo ${path}
    fi
done < "${full_species}" > /home/ninon.robin/Methode_Generic/data/sort_species_ninon.txt