#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Input / output
# ============================================================

QUERY="results_summary/high_priority_panel_external_contigs.fa"

OUTDIR="results/broad_blast"
LOGDIR="logs"

SUBMIT_RESPONSE="${OUTDIR}/ERR14789190_blast_submit.txt"
RESULT="${OUTDIR}/ERR14789190_high_priority_core_nt_megablast.csv"

API="https://blast.ncbi.nlm.nih.gov/Blast.cgi"

mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found or empty:"
    echo "$QUERY"
    exit 1
fi


# ============================================================
# Basic information
# ============================================================

NSEQ=$(grep -c "^>" "$QUERY")

echo "Query file: $QUERY"
echo "Query sequences: $NSEQ"
echo "Database: core_nt"
echo "Program: blastn / megablast"
echo


# ============================================================
# Submit one BLAST job containing all queries
# ============================================================

echo "Submitting BLAST job to NCBI..."

curl -sS \
    --retry 3 \
    --retry-delay 5 \
    --connect-timeout 30 \
    --data-urlencode "CMD=Put" \
    --data-urlencode "PROGRAM=blastn" \
    --data-urlencode "MEGABLAST=on" \
    --data-urlencode "DATABASE=core_nt" \
    --data-urlencode "EXPECT=1e-10" \
    --data-urlencode "HITLIST_SIZE=20" \
    --data-urlencode "FILTER=mL" \
    --data-urlencode "TOOL=wwshotgun_viromics" \
    --data-urlencode "QUERY@${QUERY}" \
    "$API" \
    > "$SUBMIT_RESPONSE"


# ============================================================
# Extract RID and estimated runtime
# ============================================================

RID=$(
    awk -F'=' '
    /RID =/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
    ' "$SUBMIT_RESPONSE"
)

RTOE=$(
    awk -F'=' '
    /RTOE =/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
    ' "$SUBMIT_RESPONSE"
)

if [[ -z "$RID" ]]; then
    echo "ERROR: NCBI did not return an RID."
    echo
    echo "Submission response:"
    cat "$SUBMIT_RESPONSE"
    exit 1
fi

echo "RID: $RID"

if [[ -n "$RTOE" ]]; then
    echo "NCBI estimated runtime: ${RTOE} seconds"
else
    RTOE=60
    echo "No RTOE returned; using initial 60-second wait."
fi


# ============================================================
# Initial wait
#
# Do not poll NCBI too frequently.
# ============================================================

if (( RTOE < 60 )); then
    WAIT=60
else
    WAIT=$RTOE
fi

echo "Waiting ${WAIT} seconds before first retrieval..."
sleep "$WAIT"


# ============================================================
# Poll once per minute until result is ready
# ============================================================

TMP="${OUTDIR}/.${RID}.tmp"

while true; do

    echo "Checking RID $RID ..."

    curl -sS \
        --retry 3 \
        --retry-delay 5 \
        --connect-timeout 30 \
        --data-urlencode "CMD=Get" \
        --data-urlencode "RID=${RID}" \
        --data-urlencode "FORMAT_TYPE=CSV" \
        --data-urlencode "ALIGNMENT_VIEW=Tabular" \
        --data-urlencode "DESCRIPTIONS=20" \
        --data-urlencode "ALIGNMENTS=20" \
        "$API" \
        > "$TMP"

    if grep -q "Status=WAITING" "$TMP"; then
        echo "Status: WAITING"
        sleep 60
        continue
    fi

    if grep -q "Status=UNKNOWN" "$TMP"; then
        echo "ERROR: RID is unknown or expired."
        cat "$TMP"
        exit 1
    fi

    if grep -q "Status=FAILED" "$TMP"; then
        echo "ERROR: NCBI BLAST job failed."
        cat "$TMP"
        exit 1
    fi

    # If none of the waiting/error states appears,
    # treat the returned content as the completed result.
    mv "$TMP" "$RESULT"

    echo
    echo "BLAST completed."
    echo "RID: $RID"
    echo "Result:"
    echo "$RESULT"

    break

done
