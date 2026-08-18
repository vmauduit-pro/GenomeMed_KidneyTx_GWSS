#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --export=ALL
#SBATCH --mem=50G
#SBATCH --array=22-22
#SBATCH -e logs/errlog_%x_%A_%a.out
#SBATCH -o logs/log_%x_%A_%a.out

show_help() {
  cat <<EOF
Usage:
  sbatch 1.compute_MM.sh <vcf_dir> <pairs> <outdir>
  sbatch 1.compute_MM.sh --vcf /path/to/all_chr.vcf.gz --pairs pairs.txt --outdir /path/to/out
  sbatch 1.compute_MM.sh --vcf-dir /path/to/vcfs --pairs pairs.txt --outdir /path/to/out

Positional arguments (deprecated, still supported):
  vcf_dir   Directory containing chr1.vcf.gz, chr2.vcf.gz, ... chr22.vcf.gz
  pairs     Path to pairs file (two cols: recipient donor; no header). IDs must match VCF samples.
            Supports duplicates in col1 or col2 (same recipient/donor in multiple pairs).
  outdir    Output directory (creates outdir/MM/ and outdir/PAIRS/)

Named options:
  --vcf PATH       Single VCF containing all chromosomes (preferred)
  --vcf-dir PATH   Directory containing per-chromosome chr*.vcf.gz files (legacy)
  --pairs PATH     Pairs file (recipient donor, no header) (required)
  --outdir PATH    Output directory (required)
  -h, --help       Show this help and exit
EOF
}

vcf=""
vcf_dir=""
pairs=""
outdir=""

if [[ $# -eq 0 ]]; then
  echo "Error: no arguments provided." >&2
  show_help
  exit 1
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  show_help
  exit 0
fi

if [[ "$1" != -* ]]; then
  vcf_dir=$1
  pairs=$2
  outdir=$3
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --vcf)
        vcf="$2"
        shift 2
        ;;
      --vcf-dir)
        vcf_dir="$2"
        shift 2
        ;;
      --pairs)
        pairs="$2"
        shift 2
        ;;
      --outdir)
        outdir="$2"
        shift 2
        ;;
      -h|--help)
        show_help
        exit 0
        ;;
      *)
        echo "Unknown option: $1" >&2
        show_help
        exit 1
        ;;
    esac
  done
fi

if [[ -z "$pairs" || -z "$outdir" || ( -z "$vcf" && -z "$vcf_dir" ) ]]; then
  echo "Error: one of --vcf or --vcf-dir plus --pairs and --outdir are required." >&2
  show_help
  exit 1
fi

CHR=$SLURM_ARRAY_TASK_ID

OUTDIR="${outdir%/}"

mkdir -p logs
mkdir -p "$OUTDIR/MM"
mkdir -p "$OUTDIR/PAIRS"

cd "$OUTDIR"

TMPDIR=$(mktemp -d tempdir.XXXXX)

if [[ -n "$vcf" ]]; then
  FULL_VCF="$vcf"
else
  VCF_DIR="${vcf_dir%/}"
fi

# Parse and filter pairs (exclude header-like rows, NA)
awk '{print $1,$2}' "$pairs"  | grep -v "^$" > $TMPDIR/pairs_raw.txt

# Sample list (from full VCF if provided, otherwise from chr22 file in VCF_DIR)
if [[ -n "$FULL_VCF" ]]; then
  bcftools query -l "$FULL_VCF" > $TMPDIR/all_samples.txt
else
  bcftools query -l "$VCF_DIR/chr22.vcf.gz" > $TMPDIR/all_samples.txt
fi

# Remove pairs where recipient or donor not in VCF
grep -wFvf $TMPDIR/all_samples.txt <(awk '{print $1}' $TMPDIR/pairs_raw.txt) > $TMPDIR/recipients_rm.txt
grep -wFvf $TMPDIR/all_samples.txt <(awk '{print $2}' $TMPDIR/pairs_raw.txt) > $TMPDIR/donors_rm.txt
grep -wFvf <(cat $TMPDIR/recipients_rm.txt $TMPDIR/donors_rm.txt) $TMPDIR/pairs_raw.txt > $TMPDIR/pairs.txt

# Unique donors and recipients (first-appearance order) for extraction
awk '{print $2}' $TMPDIR/pairs.txt | awk '!seen[$1]++{print}' > $TMPDIR/unique_donors.txt
awk '{print $1}' $TMPDIR/pairs.txt | awk '!seen[$1]++{print}' > $TMPDIR/unique_recipients.txt

# Build column index mapping: for each pair, which column in unique_donors/unique_recipients (1-based)
# Output: don_col_map (one index per line), rec_col_map (one per line)
awk 'NR==FNR{rec[$1]=NR; next} {print rec[$1]}' $TMPDIR/unique_recipients.txt <(awk '{print $1}' $TMPDIR/pairs.txt) > $TMPDIR/rec_col_map.txt
awk 'NR==FNR{don[$1]=NR; next} {print don[$2]}' $TMPDIR/unique_donors.txt <(awk '{print $1,$2}' $TMPDIR/pairs.txt) > $TMPDIR/don_col_map.txt

# Extract unique donors and recipients (no reheader), from full VCF or per-chromosome files
if [[ -n "$FULL_VCF" ]]; then
  bcftools view -r "$CHR" -S $TMPDIR/unique_donors.txt -Oz -o $TMPDIR/Dtmp.chr${CHR}.vcf.gz "$FULL_VCF"
else
  bcftools view -S $TMPDIR/unique_donors.txt -Oz -o $TMPDIR/Dtmp.chr${CHR}.vcf.gz "$VCF_DIR/chr${CHR}.vcf.gz"
fi
bcftools index -f $TMPDIR/Dtmp.chr${CHR}.vcf.gz

if [[ -n "$FULL_VCF" ]]; then
  bcftools view -r "$CHR" -S $TMPDIR/unique_recipients.txt -Oz -o $TMPDIR/Rtmp.chr${CHR}.vcf.gz "$FULL_VCF"
else
  bcftools view -S $TMPDIR/unique_recipients.txt -Oz -o $TMPDIR/Rtmp.chr${CHR}.vcf.gz "$VCF_DIR/chr${CHR}.vcf.gz"
fi
bcftools index -f $TMPDIR/Rtmp.chr${CHR}.vcf.gz

# Dosage matrices (columns = unique donors / unique recipients)
bcftools query -f "[%DS\t]\n" $TMPDIR/Dtmp.chr${CHR}.vcf.gz > $TMPDIR/dosage.donors.chr${CHR}.txt
bcftools query -f "[%DS\t]\n" $TMPDIR/Rtmp.chr${CHR}.vcf.gz > $TMPDIR/dosage.recipients.chr${CHR}.txt

# Build pair-ordered dosage: for each variant, output [don1, rec1, don2, rec2, ...] by column indexing
# Uses Python for clarity; works with or without duplicates in pairs
python3 - "$TMPDIR" "$CHR" << 'PYTHON_SCRIPT'
import sys
tmpdir = sys.argv[1]
chrn = sys.argv[2]

with open(f"{tmpdir}/don_col_map.txt") as f:
    don_col = [int(line.strip()) for line in f]
with open(f"{tmpdir}/rec_col_map.txt") as f:
    rec_col = [int(line.strip()) for line in f]

with open(f"{tmpdir}/dosage.donors.chr{chrn}.txt") as fd, \
     open(f"{tmpdir}/dosage.recipients.chr{chrn}.txt") as fr, \
     open(f"{tmpdir}/dosage.chr{chrn}.txt", "w") as out:
    for don_line, rec_line in zip(fd, fr):
        don_vals = don_line.rstrip().split("\t")
        rec_vals = rec_line.rstrip().split("\t")
        parts = []
        for i in range(len(don_col)):
            parts.append(don_vals[don_col[i]-1])
            parts.append(rec_vals[rec_col[i]-1])
        out.write("\t".join(parts) + "\n")
PYTHON_SCRIPT

# Compute mismatches (x=donor, y=recipient). Columns are [don1, rec1, don2, rec2, ...]
awk '
function abs(x) { return x < 0 ? -x : x }
{
    npairs = NF/2
    for (j = 1; j <= npairs; j++) {
        x = int($(2*j-1) + 0.5)   # donor j
        y = int($(2*j) + 0.5)     # recipient j
        result = abs(y - 1) * abs(x - y)
        printf result"\t"
    }
    print ""
}
' $TMPDIR/dosage.chr${CHR}.txt | sed "s/\t$//g" > $TMPDIR/MM.chr${CHR}.txt

# Reconstruct VCF: header with sample names recipient_donor (pair order), body from MM matrix
# Use first available VCF as template for #CHROM..FORMAT; sample names from pairs
sample_names=$(awk '{print $1"_"$2}' $TMPDIR/pairs.txt | paste -sd$'\t')
# Build minimal header
bcftools view -h $TMPDIR/Dtmp.chr${CHR}.vcf.gz | head -n -1 > $TMPDIR/MM.header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_names" >> $TMPDIR/MM.header.txt

# Format MM values as DS (dosage) for plink2 compatibility
bcftools view -H $TMPDIR/Dtmp.chr${CHR}.vcf.gz | cut -f1-9 | sed "s/GT:DS:HDS:GP/DS/g" > $TMPDIR/MM.meta.txt
# Body: meta + MM values as DS (one value per sample, semicolon not needed for single value)
awk 'NR==FNR{meta[NR]=$0; next} {
    line = meta[FNR]
    for (i=1; i<=NF; i++) line = line "\t" $i
    print line
}' $TMPDIR/MM.meta.txt $TMPDIR/MM.chr${CHR}.txt > $TMPDIR/MM.body.chr${CHR}.txt

# Full VCF (header without duplicate FORMAT IDs + body)
(head -n -1 $TMPDIR/MM.header.txt | grep -v "ID=GP" | grep -v "ID=GT" | grep -v "ID=HDS"; tail -1 $TMPDIR/MM.header.txt; cat $TMPDIR/MM.body.chr${CHR}.txt) > $TMPDIR/chr${CHR}.vcf

bcftools sort -T . $TMPDIR/chr${CHR}.vcf -Oz -o MM/chr${CHR}.vcf.gz
bcftools index -f MM/chr${CHR}.vcf.gz

plink2 --vcf $TMPDIR/chr${CHR}.vcf dosage=DS --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 --make-pgen --out MM/chr${CHR}

if [ "$CHR" -eq 22 ]; then
    cp $TMPDIR/pairs.txt PAIRS/full_pairs.txt
fi

rm -rf $TMPDIR
