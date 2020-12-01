# Download the GRCh38 annotation file

mkdir -p ../GRCh38
cd  ../GRCh38
wget -nc -P ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz

# If need be, unzip
gunzip gencode.v26.annotation.gtf.gz
