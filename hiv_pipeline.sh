#!/usr/bin/env bash
# ----------------------------------------------------------------------
# HIV Drug Resistance Pipeline v0.1
# Author: Ronan Doyle, Lead Clinical Bioinformatician, Synnovis
# ----------------------------------------------------------------------
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="/data/HIV_pipeline_results"

# ============================================================== 
# ARG PARSE
# ============================================================== 

THREADS=1
SAMPLESHEET=""
RUN_DIR=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --samplesheet) SAMPLESHEET="$2"; shift; shift ;;
        --threads) THREADS="$2"; shift; shift ;;
        --rundir) RUN_DIR="$2"; shift; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$SAMPLESHEET" || -z "$THREADS" || -z "$RUN_DIR" ]]; then
    echo "❌ ERROR: Must specify --samplesheet, --threads, --rundir"
    exit 1
fi

RUN_DIR=$(realpath "$RUN_DIR")
[[ -d "$RUN_DIR" ]] || { echo "Run directory not found"; exit 1; }

SAMPLESHEET=$(realpath "$SAMPLESHEET")
[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found"; exit 1; }

RUN_NAME=$(basename "$RUN_DIR")

# ============================================================== 
# PIPELINE LOG DIR
# ============================================================== 

mkdir -p "${RESULTS_DIR}/pipeline_logs"
RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

touch "${RESULTS_DIR}/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"
MASTER_SUMMARY="${RESULTS_DIR}/pipeline_logs/pipeline_summary_${RUN_NAME}_${RUN_TIMESTAMP}.txt"

echo -e "Sample\tStatus\tFailed_Step\tTotal_Time(s)" > "$MASTER_SUMMARY"

# ============================================================== 
# READ SAMPLESHEET
# ============================================================== 

SAMPLES=(); BARCODES=()
while IFS=, read -r barcode sample || [[ -n "$barcode" ]]; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    SAMPLES+=("$sample")
    BARCODES+=("$barcode")
done < "$SAMPLESHEET"

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 

process_sample() {
    local idx="$1"
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"

    sample_dir="${RESULTS_DIR}/${sample}"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1
    trap 'exec 1>&3 2>&4; exec 3>&- 4>&-' RETURN

    log "========== Processing sample: $sample (Barcode $barcode) =========="
    start_total=$(date +%s)
    status="OK"; failed_step="-"

    # ---------------- STEP 1: Concatenate ----------------
    step="concat"; step_start=$(date +%s)
    log "▶ Step $step"
    src_pattern="${RUN_DIR}/**/fastq_pass/barcode${barcode}/*.fastq.gz"
    shopt -s globstar nullglob
    files=( $src_pattern )
    if (( ${#files[@]} == 0 )); then
        log "⚠️ No FASTQ found for barcode${barcode}"
        status="FAILED"; failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        return
    fi
    log "Concatenating ${#files[@]} fastq files → ${sample_dir}/${sample}_raw.fastq.gz"
    if ! cat "${files[@]}" > "${sample_dir}/${sample}_raw.fastq.gz"; then
       log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 2: Porechop & NanoStat ----------------
    step="porechop_nanostat"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/qc
    fq="${sample_dir}/${sample}_raw.fastq.gz"; out_trim="${sample_dir}/${sample}.fastq.gz"
    if ! porechop -i "$fq" -o "$out_trim" -t "$THREADS" --no_split; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS"
    rm -vf "$fq"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Fastq2CodFreq ----------------
    step="fastq2codfreq"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! fastq2codfreq -p "minimap2" -r "${SCRIPT_DIR}/profiles/HIV1.json" --workers "$THREADS" --no-autopairing "${sample_dir}/"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 4: HIVDB SeqReads analysis ----------------
    step="seqreads"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! sierrapy seqreads -D "50" -m "0.02" "${sample_dir}/${sample}.codfreq"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: Produce Report ----------------
    step="html_report"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! python "${SCRIPT_DIR}/htmlview/HTMLMaker.py" -o "${sample_dir}/${sample}.html" "${sample_dir}/${sample}.report.json"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="ANALYSIS COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample analysis completed in ${total_time}s"
}

for idx in "${!SAMPLES[@]}"; do
    process_sample "$idx"
done

echo "Pipeline complete. Run log saved in $MASTER_SUMMARY"