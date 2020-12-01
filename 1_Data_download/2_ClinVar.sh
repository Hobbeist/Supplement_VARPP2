# Download the ClinVar annotation of known pathogenic variants

wget -nc -P ../data https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

#Unzip
unzip clinvar.vcf.gz

# Extract level 5 clinical significant genes: 'Pathogenic'
cat ../data/clinvar.vcf | grep 'CLNSIG' | grep 'Pathogenic' | awk {'print $1 , $2 , $3 , $4'} > ../data/ClinVar_pathogenic
sed -e'1i\CHROM\tPOS\tREF\tALT' ../data/ClinVar_pathogenic > ../data/ClinVar_patho
rm ../data/ClinVar_pathogenic
