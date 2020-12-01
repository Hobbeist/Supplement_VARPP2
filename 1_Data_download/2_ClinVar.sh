# Download the ClinVar annotation of known pathogenic variants

mkdir -p ../ClinVar
wget -nc -P ../ClinVar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

#Unzip
unzip clinvar.vcf.gz

# Extract level 5 clinical significant genes: 'Pathogenic'
cat ../ClinVar/clinvar.vcf | grep 'CLNSIG' | grep 'Pathogenic' | awk {'print $1 , $2 , $3 , $4'} > ../ClinVar/ClinVar_pathogenic
sed -e'1i\CHROM\tPOS\tREF\tALT' ../ClinVar/ClinVar_pathogenic > ../ClinVar/ClinVar_patho
rm ../ClinVar/ClinVar_pathogenic
