#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:59:00

SAMPLE="Lakota"

module load samtools
module load bcftools
module load picard
module load bwa

mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/bam_dir
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/vcf_files
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/vcf_files/scaffolds
mkdir /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE/bam_dir/bam_scaffs

REF="/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_Ovbor_1.2.fna" #CHANGE IF SAMPLE IS A MD
RAW="/home/devan/projects/rrg-shaferab/DEER_FASTQ/SPERM/HiC/AS-11671/$SAMPLE"
RAW2="/home/devan/projects/rrg-shaferab/DEER_FASTQ/SPERM/HiC2/AS-11671/$SAMPLE"
MAIN="/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC/$SAMPLE"
BAMDIR="$MAIN/bam_dir"
VCFDIR="$MAIN/vcf_files"

# Align to a reference genome using bwa mem and importantly the -5SP -T0 options
echo "beginning mapping for "$RAW/$SAMPLE"_R1.fastq.gz and "$RAW/$SAMPLE"_R2.fastq.gz"
bwa mem -5SP -T0 "$REF" "$RAW/${SAMPLE}_R1.fastq.gz" "$RAW/${SAMPLE}_R2.fastq.gz" | \
samtools view -Sb | samtools sort -o "$BAMDIR/$SAMPLE.set1.sorted.bam"
echo "mapping done for run 1"

echo "beginning mapping for "$RAW2/$SAMPLE"_R1.fastq.gz and "$RAW2/$SAMPLE"_R2.fastq.gz"
bwa mem -5SP -T0 "$REF" "$RAW2/${SAMPLE}_R1.fastq.gz" "$RAW2/${SAMPLE}_R2.fastq.gz" | \
samtools view -Sb - | samtools sort - -o "$BAMDIR/$SAMPLE.set2.sorted.bam"
echo "mapping done for run 2"

# Merge and mark duplicates with picard
echo "merging and marking duplicates"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I="$BAMDIR/$SAMPLE.set1.sorted.bam" \
      I="$BAMDIR/$SAMPLE.set2.sorted.bam" \
      O="$BAMDIR/$SAMPLE.merged.mkdup.bam" \
      M="$BAMDIR/$SAMPLE.metrics.txt" \
      CREATE_INDEX=true

# Call variants directly using the Hi-C reads with bcftools.
# (Stage 1)
bcftools mpileup --count-orphans -Ou -f $REF "${BAMDIR}/$SAMPLE.merged.mkdup.bam" -o "${VCFDIR}/$SAMPLE.mpileup.bcf"
# (Stage 2)
bcftools call --threads 8 -mv -Oz "${VCFDIR}/$SAMPLE.mpileup.bcf" -o "${VCFDIR}/$SAMPLE.unfiltered.vcf.gz"
bcftools index "${VCFDIR}/$SAMPLE.unfiltered.vcf.gz"
# Define input VCF file
INPUT_VCF="${VCFDIR}/$SAMPLE.unfiltered.vcf.gz" 

# Loop through each scaffold name in the list file
while IFS= read -r SCAFFOLD_ID; do
  
  # Set the output filename
  OUTPUT_VCF="${VCFDIR}/scaffolds/${SCAFFOLD_ID}.vcf"
  
  echo "Processing $SCAFFOLD_ID..."
  
  # Use bcftools view to filter by region (-r)
  bcftools view -r "$SCAFFOLD_ID" \
    -o "$OUTPUT_VCF" \
    -O v \
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
  samtools view -b -@ 4 -o "$OUTPUT_BAM" "${BAMDIR}/$SAMPLE.merged.mkdup.bam" "$SCAFFOLD_ID"
  samtools index "$OUTPUT_BAM"
    
done < WTD_L90_scaffs.list
