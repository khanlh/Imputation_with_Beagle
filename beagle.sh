#!/bin/bash

set -euo pipefail  # Exit if any command fails, if any variable is unset, or in pipeline failure

# ===== INPUT =====
REF="/data/for_others/genome_index/genome.fa"  # Reference genome FASTA
HIGHCOV="/data/KhanLe_Data/workshop_lobi/ws3/bam_wo_trim/bcftools/vcf_output/merged.vcf.gz"  # High-coverage multi-sample VCF
LOWCOV_LIST="/data/KhanLe_Data/vannamei/lowpass/0.5X/vcf/*.vcf.gz"  # List of low-coverage sample VCF files

# ===== STEP 1: Normalize high-coverage reference panel =====
echo "[Step 1] Normalizing high-coverage panel..."

source /data/conda_envs/bcftools/envs/bin/activate
export BCFTOOLS_PLUGINS=/data/conda_envs/bcftools/envs/libexec/bcftools

mkdir -p /data/KhanLe_Data/vannamei/lowpass/beagle

# Normalize and filter SNPs
bcftools norm -f "$REF" -m both "$HIGHCOV" -Oz -o beagle/1.highcov.norm.vcf.gz
bcftools view -v snps -m2 -M2 beagle/1.highcov.norm.vcf.gz -Oz -o beagle/2.highcov.snps.vcf.gz

# Add MAF and missing rate tags, then filter
bcftools +fill-tags beagle/2.highcov.snps.vcf.gz -- -t MAF,F_MISSING | \
    bcftools view -i 'MAF>0.01 && F_MISSING<0.05' -Oz -o beagle/3.highcov.snps.qc.vcf.gz

# Index the filtered VCF
bcftools index beagle/3.highcov.snps.qc.vcf.gz

# Remove missing genotypes
bcftools view -g ^miss beagle/3.highcov.snps.qc.vcf.gz -Oz -o beagle/3.highcov.nomiss.vcf.gz

# Convert unphased genotypes to phased format
bcftools view beagle/3.highcov.nomiss.vcf.gz | \
    sed 's/0\/0/0|0/g; s/0\/1/0|1/g; s/1\/0/1|0/g; s/1\/1/1|1/g' | \
    bgzip -c > beagle/3.highcov.phased.vcf.gz

# Index the phased VCF
tabix -p vcf beagle/3.highcov.phased.vcf.gz

# ===== STEP 2: Normalize each low-coverage sample =====
echo "[Step 2] Normalizing low-coverage samples..."

for f in $LOWCOV_LIST; do
    base=$(basename "$f" .vcf.gz)
    echo "  - Processing sample: $base"

    bcftools norm -f "$REF" -m both "$f" -Oz -o beagle/4.${base}.norm.vcf.gz
    bcftools view -v snps -m2 -M2 beagle/4.${base}.norm.vcf.gz -Oz -o beagle/5.${base}.snps.vcf.gz
    bcftools index beagle/5.${base}.snps.vcf.gz
done

# Merge all normalized low-coverage SNP VCFs
echo "[Step 2b] Merging low-coverage samples..."
bcftools merge beagle/5.*.snps.vcf.gz -Oz -o beagle/6.lowcov.merged.snps.vcf.gz
bcftools index beagle/6.lowcov.merged.snps.vcf.gz

# ===== STEP 3: Imputation =====
echo "[Step 3] Running Beagle for imputation..."

java -Xmx32g -jar /usr/Beagle/beagle.01Mar24.d36.jar \
    gt=beagle/6.lowcov.merged.snps.vcf.gz \
    ref=beagle/3.highcov.phased.vcf.gz \
    out=beagle/7.lowcov.imputed.vcf.gz \
    nthreads=32

# ===== STEP 4: Accuracy evaluation =====
echo "[Step 4] Evaluating imputation accuracy..."

for s in $(bcftools query -l beagle/7.lowcov.imputed.vcf.gz); do
    echo "  -> Evaluating sample: $s"

    # Extract ground truth SNPs from high-coverage VCF
    bcftools view -s "$s" beagle/3.highcov.snps.qc.vcf.gz -Oz -o beagle/8.${s}.highcov.vcf.gz
    bcftools index beagle/8.${s}.highcov.vcf.gz

    # Extract imputed SNPs
    bcftools view -s "$s" beagle/7.lowcov.imputed.vcf.gz -Oz -o beagle/9.${s}.imputed.vcf.gz
    bcftools index beagle/9.${s}.imputed.vcf.gz

    # Genotype concordance check
    bcftools gtcheck -g beagle/8.${s}.highcov.vcf.gz beagle/9.${s}.imputed.vcf.gz \
        > beagle/10.${s}.concordance.txt

    # Dosage R² evaluation (requires gwas-compare plugin in vcftools)
    vcftools --vcf beagle/9.${s}.imputed.vcf.gz \
        --gwas-compare beagle/8.${s}.highcov.vcf.gz \
        --out beagle/10.${s}.eval || true
done

echo "DONE!"
