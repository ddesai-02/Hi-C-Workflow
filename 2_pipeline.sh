#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem 16G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:59:00 

# =====================================
# Hi-C Recombination Pipeline
# =====================================

# --- 1. Define Paths ---
MAIN="/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/sperm/HiC"

SAMPLE_PREFIX="Lakota" # Change per individual
SCAFFOLD_LIST="$MAIN/WTD_L90_scaffs.list" # Change per species
HAPCUT2_DIR="$MAIN/HapCUT2/build"
HI_RECOMB_DIR="$MAIN/Hi-reComb/Build"


# --- 2. Define Input Directories ---
# Input files should be separated by scaffold
INPUT_BAM_DIR="$MAIN/$SAMPLE_PREFIX/bam_dir/bam_scaffs" # For HapCUT2 extractHAIRS
INPUT_VCF_DIR="$MAIN/$SAMPLE_PREFIX/vcf_files/scaffolds"  # For HapCUT2 extractHAIRS and HAPCUT2 phasing

# --- 3. Define Output Directories ---
# Create the main folder and sub-directories
OUTPUT_ROOT_DIR="$MAIN/$SAMPLE_PREFIX"
OUTPUT_FRAG_DIR="${OUTPUT_ROOT_DIR}/hapcut_fragments"       # Output from extractHAIRS
OUTPUT_HAPS_DIR="${OUTPUT_ROOT_DIR}/phased_haplotypes"      # Output from HAPCUT2
OUTPUT_SAM_DIR="${OUTPUT_ROOT_DIR}/informative_sam"         # Output from FindInfoPairs
OUTPUT_MAPS_DIR="${OUTPUT_ROOT_DIR}/recombination_maps"     # Output from RecombMap

mkdir -p "$OUTPUT_FRAG_DIR" "$OUTPUT_HAPS_DIR" "$OUTPUT_SAM_DIR" "$OUTPUT_MAPS_DIR"

echo "====================================================================="
echo "Pipeline Start: Processing scaffolds from $SCAFFOLD_LIST"
echo "Output will be in the '$OUTPUT_ROOT_DIR/' directory."
echo "====================================================================="

# --- 4. Main Loop: Iterate through each Scaffold ID ---
while IFS= read -r SCAFFOLD_ID; do
    echo -e "\n--- Starting Processing for Scaffold: **$SCAFFOLD_ID** ---"

    # Define all required file paths for the current scaffold
    BAM_IN="${INPUT_BAM_DIR}/${SCAFFOLD_ID}.bam"
    VCF_IN="${INPUT_VCF_DIR}/${SCAFFOLD_ID}.vcf"

    FRAG_OUT="${OUTPUT_FRAG_DIR}/hapcutFragments_${SCAFFOLD_ID}.txt"
    HAPS_OUT="${OUTPUT_HAPS_DIR}/${SCAFFOLD_ID}_PHASE.txt"
    SAM_OUT="${OUTPUT_SAM_DIR}/${SCAFFOLD_ID}_INFORMATIVE_READS.sam"
    MAP_OUT="${OUTPUT_MAPS_DIR}/${SCAFFOLD_ID}_RecombMap.txt"

    # =================================================================
    # STEP 1: HapCUT2 extractHAIRS (Generate Fragment File)
    # =================================================================

    "${HAPCUT2_DIR}/extractHAIRS" \
        --hic 1 \
        --VCF "$VCF_IN" \
        --bam "$BAM_IN" \
        --out "$FRAG_OUT" 2>/dev/null
    
    if [ $? -ne 0 ] || [ ! -s "$FRAG_OUT" ]; then
        echo "❌ ERROR in extractHAIRS for $SCAFFOLD_ID. Skipping scaffold."
        continue # Move to the next scaffold
    fi
    echo "   ✅ Fragments created: $FRAG_OUT"

    # =================================================================
    # STEP 2: HapCUT2 HAPCUT2 (Phase Haplotypes)
    # =================================================================
    echo "2/4. Running **HAPCUT2** (Phasing)..."
    "${HAPCUT2_DIR}/HAPCUT2" \
        --VCF "$VCF_IN" \
        --hic 1 \
        --fragments "$FRAG_OUT" \
        --output "$HAPS_OUT" 2>/dev/null

    if [ $? -ne 0 ] || [ ! -s "$HAPS_OUT" ]; then
        echo "❌ ERROR in HAPCUT2 for $SCAFFOLD_ID. Skipping scaffold."
        continue
    fi
    echo "   ✅ Phased haplotypes created: $HAPS_OUT"

    # =================================================================
    # STEP 3: Hi-reComb FindInfoPairs (Identify Informative Reads)
    # =================================================================
    echo "3/4. Running **FindInfoPairs**..."
    # samtools view -F 3332 filters out: Unmapped (4), Secondary (256), Duplicate (1024), Supplementary (2048)
    samtools view -F 3332 "$BAM_IN" "$SCAFFOLD_ID" | \
    "${HI_RECOMB_DIR}/Hi-reComb" FindInfoPairs -m 10 "$HAPS_OUT" > "$SAM_OUT"

    if [ ! -s "$SAM_OUT" ]; then
        echo "❌ ERROR in FindInfoPairs for $SCAFFOLD_ID (Output SAM is empty). Skipping recomb."
        continue
    fi
    echo "   ✅ Informative SAM created: $SAM_OUT"

    # =================================================================
    # STEP 4: Hi-reComb RecombMap (Estimate Recombination Map)
    # =================================================================
    echo "4/4. Running **RecombMap**..."
    "${HI_RECOMB_DIR}/Hi-reComb" RecombMap -d 500 -q 20 \
        "$HAPS_OUT" \
        "$SAM_OUT" \
        > "$MAP_OUT" 2>/dev/null

    if [ -s "$MAP_OUT" ]; then
        echo "   ✅ **SUCCESS**: Recombination map created: $MAP_OUT"
    else
        echo "   ⚠️ WARNING: Map file is empty or missing for $SCAFFOLD_ID."
    fi

done < "$SCAFFOLD_LIST"

echo -e "\n====================================================================="
echo "Pipeline End: All steps complete. Check '$OUTPUT_ROOT_DIR/' for results."
echo "====================================================================="
