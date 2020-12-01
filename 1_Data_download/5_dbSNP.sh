# Download the dbSNP database for all known human variants, benign and pathogenic

# Make a dbSNP directory
mkdir -p ../dbSNP
cd ../dbSNP
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
gunzip 00-All.vcf.gz
