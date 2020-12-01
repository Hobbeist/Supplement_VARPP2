# Download the dbSNP database for all known human variants, benign and pathogenic

cd ../data
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
gunzip 00-All.vcf.gz
