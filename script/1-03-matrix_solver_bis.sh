# Methode d'assemblage de la matrice si les especes ont ete traitees en n fractions :

centro='/home/ninon.robin/Methode_Generic/output/centro.tsv'

nom_fraction_1_species='/home/ninon.robin/Methode_Generic/output/species_*nom_fraction*.tsv'
#... (<- Dupliquer ici la ligne precedente n-2 fois)
nom_fraction_n_species='/home/ninon.robin/Methode_Generic/output/species_*nom_fraction*.tsv'

nom_fraction_1_share='/home/ninon.robin/Methode_Generic/output/share_*nom_fraction*.tsv'
#... (<- Dupliquer ici la ligne precedente n-2 fois)
nom_fraction_n_share='/home/ninon.robin/Methode_Generic/output/share_*nom_fraction*.tsv'

cat ${nom_fraction_1_species} ... ${nom_fraction_n_species} > /home/ninon.robin/Methode_Generic/output/matrix.tsv
paste ${centro} ${nom_fraction_1_share} ... ${nom_fraction_n_share} | sed 's/\t//g' >> /home/ninon.robin/Methode_Generic/output/matrix.tsv
