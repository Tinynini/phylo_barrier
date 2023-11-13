centro='/home/ninon.robin/Methode_Generic/output/centro.tsv'
species='/home/ninon.robin/Methode_Generic/output/species.tsv'
share='/home/ninon.robin/Methode_Generic/output/share.tsv'

cat ${species} > /home/ninon.robin/Methode_Generic/output/matrix.tsv
paste ${centro} ${share} | sed 's/\t//g' >> /home/ninon.robin/Methode_Generic/output/matrix.tsv