#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --time=4-23:59:00

SAMPLE="Lakota"

module load samblaster
module load samtools
module load bcftools

mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/bam_dir
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/vcf_files
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/vcf_files/scaffolds
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/bam_dir/bam_scaffs

REF="/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_Ovbor_1.2.fna" # CHANGE IF SAMPLE IS A MD
RAW="/home/devan/projects/rrg-shaferab/DEER_FASTQ/SPERM/HiC/AS-11671/$SAMPLE"
RAW2="/home/devan/projects/rrg-shaferab/DEER_FASTQ/SPERM/HiC2/AS-11671/$SAMPLE"
MAIN="/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE"
BAMDIR="$MAIN/bam_dir"
VCFDIR="$MAIN/vcf_files"

# Align to a reference genome using bwa mem using the -5SP -T0 options, mark duplicates with samblaster
echo "beginning mapping for "$RAW/$SAMPLE"_R1.fastq.gz and "$RAW/$SAMPLE"_R2.fastq.gz"
bwa mem -5SP -T0 "$REF" "$RAW/$SAMPLE"_R1.fastq.gz "$RAW/$SAMPLE"_R2.fastq.gz | samblaster | samtools view -S -h -b -F 2304 | samtools sort -@ 16 -n > "$BAMDIR/$SAMPLE".run1.bam
echo "mapping done for run 1"

bwa mem -5SP -T0 "$REF" "$RAW2/$SAMPLE"_R1.fastq.gz "$RAW2/$SAMPLE"_R2.fastq.gz -o "$BAMDIR/$SAMPLE".run2.tmp.bam
samblaster -i "$BAMDIR/$SAMPLE".run2.tmp.bam -o "$BAMDIR/$SAMPLE".run2.mrkdup.bam
samtools view -S -h -b -F 2304 "$BAMDIR/$SAMPLE".run2.mrkdup.bam -o "$BAMDIR/$SAMPLE".run2.tmp2.bam
samtools sort -@ 32 -n "$BAMDIR/$SAMPLE".run2.tmp2.bam -o "$BAMDIR/$SAMPLE".run2.bam
echo "mapping done for run 2"

samtools merge -@16 -n -o "${BAMDIR}/${SAMPLE}.merged.bam" "${BAMDIR}/${SAMPLE}.run1.bam" "${BAMDIR}/${SAMPLE}.run2.bam"
echo "merging done"

samtools sort -@ 16 -o "${BAMDIR}/${SAMPLE}.merged.sorted.bam" "${BAMDIR}/${SAMPLE}.merged.bam"
samtools index "${BAMDIR}/${SAMPLE}.merged.sorted.bam"

# Call variants directly using the Hi-C reads with bcftools.
# (Stage 1)
bcftools mpileup --count-orphans -Ou -f $REF "${BAMDIR}/${SAMPLE}.merged.sorted.bam" -o "${VCFDIR}/$SAMPLE.mpileup.bcf"
# (Stage 2)
bcftools call --threads 8 -mv -Oz "${VCFDIR}/$SAMPLE.mpileup.bcf" -o "${VCFDIR}/$SAMPLE.unfiltered.vcf.gz"
# (Stage 3) - Index
tabix "${VCFDIR}/$SAMPLE.unfiltered.vcf.gz"

# Define your input VCF file
INPUT_VCF="${VCFDIR}/$SAMPLE.unfiltered.vcf.gz" 

# Loop through each scaffold name in the list file
while IFS= read -r SCAFFOLD_ID; do
  
  # Set the output filename
  OUTPUT_VCF="${VCFDIR}/scaffolds/${SCAFFOLD_ID}.vcf.gz"
  
  echo "Processing $SCAFFOLD_ID..."
  
  # Use bcftools view to filter by region (-r)
  bcftools view -r "$SCAFFOLD_ID" \
    -o "$OUTPUT_VCF" \
    -O z \
    "$INPUT_VCF"
    
done < WTD_L90_scaffs.list

# Loop through each scaffold name in the list file
while IFS= read -r SCAFFOLD_ID; do
  
  # Set the output filename
  OUTPUT_BAM="${BAMDIR}/bam_scaffs/${SCAFFOLD_ID}.bam"
  
  echo "Processing $SCAFFOLD_ID..."
  
  # Use samtools view to filter by region (Scaffold ID)
  # -b: Output in BAM format
  # -@ 4: Use 4 threads (adjust as needed for faster processing)
  samtools view -b -@ 4 -o "$OUTPUT_BAM" "${BAMDIR}/$SAMPLE.merged.sorted.bam" "$SCAFFOLD_ID"
    
done < WTD_L90_scaffs.list
