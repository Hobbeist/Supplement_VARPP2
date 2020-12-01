# We need to create bedfiles out of all the data we have downloaded. This is done so we can use bedtools, the
# exteremly efficient tool to merge and intersect files based on gene coordinate

cd ../data

# CLinVar
cat ClinVar_patho | sed '1d' | awk {'print $1="chr"$1, $2, $2'} | tr ' ' '\t' > ClinVar.bed

# CADD
cat whole_genome_SNVs.tsv | sed '1,2d' > CADD.bed

# Here we include the 5th column, as this is the CADD score. We need to merge this into the
# intersect data
# The $2-1 is necessary as for bedtools we need distinct positions. Even if it is just a single base, we need the start and end. and that would 
# be start = start -1
cat CADD.bed | awk {'print $1="chr"$1, $2-1, $2, $5'} | tr ' ' '\t'  > CADD_final.bed


#dbSNP

# Extract Chromosome, position, Ref and Alt
../bcftools/bcftools query -f '%CHROM %POS %REF %ALT\n' 00-All.vcf > dbSNP.bed

# Create the bedfile with the necessary fromat.
cat dbSNP.bed | awk {'print $1="chr"$1, $2-1, $2'} | tr ' ' '\t' > dbSNP_final.bed

# remove redundant file
rm dbSNP.bed

# Sort the bedfile for intersect
bedtools sort -i dbSNP.bed > dbSNP_sorted.bed

