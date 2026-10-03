#!/usr/bin/env bash

set -u -o pipefail

IDS="results_summary/protein_level_candidates.txt"
SOURCE_FA="results_summary/protein_level_candidates.fa"

BATCHDIR="results_summary/blastx_batches"
OUTDIR="results/protein_blast"
LOGDIR="logs"

API="https://blast.ncbi.nlm.nih.gov/Blast.cgi"

COMBINED="${OUTDIR}/ERR14789190_blastx_clusterednr.csv"
SUMMARY="${OUTDIR}/ERR14789190_blastx_jobs.tsv"

BATCH_SIZE=3
POLL_INTERVAL=60
MAX_POLLS=180

mkdir -p "$BATCHDIR" "$OUTDIR" "$LOGDIR"

if [[ ! -s "$IDS" || ! -s "$SOURCE_FA" ]]; then
    echo "ERROR: Missing input files."
    exit 1
fi

echo "Protein-level candidate contigs: $(wc -l < "$IDS")"
echo "Program: BLASTx"
echo "Database: ClusteredNR (nr_cluster_seq)"
echo

# ------------------------------------------------------------
# Prepare batches
# ------------------------------------------------------------

rm -f "${BATCHDIR}"/batch_*

split \
    -l "$BATCH_SIZE" \
    -d \
    -a 2 \
    "$IDS" \
    "${BATCHDIR}/batch_"

for ids_file in "${BATCHDIR}"/batch_[0-9][0-9]; do

    batch=$(basename "$ids_file")

    awk '
    NR==FNR {
        keep[$1]=1
        next
    }
    /^>/ {
        id=substr($1,2)
        printseq=(id in keep)
    }
    printseq
    ' \
    "$ids_file" \
    "$SOURCE_FA" \
    > "${BATCHDIR}/${batch}.fa"

done

echo "Prepared:"
for fa in "${BATCHDIR}"/batch_*.fa; do
    echo "$(basename "$fa"): $(grep -c "^>" "$fa") sequences"
done
echo

: > "$COMBINED"
printf "batch\tRID\tn_queries\tstatus\n" > "$SUMMARY"


# ------------------------------------------------------------
# Sequential BLASTx
# ------------------------------------------------------------

for QUERY in "${BATCHDIR}"/batch_*.fa; do

    batch=$(basename "$QUERY" .fa)
    nq=$(grep -c "^>" "$QUERY")

    SUBMIT="${OUTDIR}/${batch}_blastx_submit.txt"
    STATUS_FILE="${OUTDIR}/.${batch}_blastx_status.tmp"
    RESULT_TMP="${OUTDIR}/.${batch}_blastx_result.tmp"
    RESULT="${BATCHDIR}/${batch}.blastx.csv"

    echo "============================================================"
    echo "Starting $batch ($nq queries)"
    echo "============================================================"

    if ! curl -sS \
        --retry 3 \
        --retry-delay 10 \
        --connect-timeout 60 \
        --max-time 180 \
        --data-urlencode "CMD=Put" \
        --data-urlencode "PROGRAM=blastx" \
        --data-urlencode "DATABASE=nr_cluster_seq" \
        --data-urlencode "EXPECT=1e-5" \
        --data-urlencode "HITLIST_SIZE=20" \
        --data-urlencode "TOOL=wwshotgun_viromics" \
        --data-urlencode "QUERY@${QUERY}" \
        "$API" \
        > "$SUBMIT"
    then
        echo "Submission failed."
        printf "%s\tNA\t%s\tSUBMISSION_FAILED\n" \
            "$batch" "$nq" >> "$SUMMARY"
        continue
    fi

    RID=$(
        awk -F'=' '
        /RID =/ {
            gsub(/[[:space:]]/, "", $2)
            print $2
            exit
        }' "$SUBMIT"
    )

    RTOE=$(
        awk -F'=' '
        /RTOE =/ {
            gsub(/[[:space:]]/, "", $2)
            print $2
            exit
        }' "$SUBMIT"
    )

    if [[ -z "$RID" ]]; then
        echo "No RID returned."
        printf "%s\tNA\t%s\tNO_RID\n" \
            "$batch" "$nq" >> "$SUMMARY"
        continue
    fi

    echo "RID: $RID"

    [[ -z "$RTOE" ]] && RTOE=60
    (( RTOE < 60 )) && WAIT=60 || WAIT=$RTOE

    echo "Initial wait: ${WAIT}s"
    sleep "$WAIT"

    success=0

    for ((poll=1; poll<=MAX_POLLS; poll++)); do

        echo "$(date) Checking RID $RID ..."

        if ! curl -sS \
            --retry 3 \
            --retry-delay 10 \
            --connect-timeout 60 \
            --max-time 120 \
            --data-urlencode "CMD=Get" \
            --data-urlencode "RID=${RID}" \
            --data-urlencode "FORMAT_OBJECT=SearchInfo" \
            "$API" \
            > "$STATUS_FILE"
        then
            echo "Network error; retrying..."
            sleep "$POLL_INTERVAL"
            continue
        fi

        if grep -qiE \
            'CPU usage limit was exceeded|SIGXCPU|An error has occurred' \
            "$STATUS_FILE"
        then
            echo "NCBI server error."
            printf "%s\t%s\t%s\tSERVER_ERROR\n" \
                "$batch" "$RID" "$nq" >> "$SUMMARY"
            break
        fi

        STATUS=$(
            grep -oE 'Status=(WAITING|READY|FAILED|UNKNOWN)' \
            "$STATUS_FILE" \
            | head -1 \
            | cut -d= -f2
        )

        if [[ "$STATUS" == "WAITING" ]]; then
            echo "Status: WAITING"
            sleep "$POLL_INTERVAL"
            continue
        fi

        if [[ "$STATUS" == "FAILED" || "$STATUS" == "UNKNOWN" ]]; then
            echo "Status: $STATUS"
            printf "%s\t%s\t%s\t%s\n" \
                "$batch" "$RID" "$nq" "$STATUS" >> "$SUMMARY"
            break
        fi

        if [[ "$STATUS" == "READY" ]]; then

            echo "Status: READY"

            if ! curl -sS \
                --retry 3 \
                --retry-delay 10 \
                --connect-timeout 60 \
                --max-time 180 \
                --data-urlencode "CMD=Get" \
                --data-urlencode "RID=${RID}" \
                --data-urlencode "FORMAT_TYPE=CSV" \
                --data-urlencode "ALIGNMENT_VIEW=Tabular" \
                --data-urlencode "DESCRIPTIONS=20" \
                --data-urlencode "ALIGNMENTS=20" \
                "$API" \
                > "$RESULT_TMP"
            then
                echo "Download error; retrying..."
                sleep "$POLL_INTERVAL"
                continue
            fi

            if grep -qiE \
                'CPU usage limit was exceeded|SIGXCPU|An error has occurred' \
                "$RESULT_TMP"
            then
                echo "Result contains NCBI error."
                printf "%s\t%s\t%s\tRESULT_ERROR\n" \
                    "$batch" "$RID" "$nq" >> "$SUMMARY"
                break
            fi

            mv "$RESULT_TMP" "$RESULT"

            # Add only non-empty data lines
            sed '/^[[:space:]]*$/d' "$RESULT" >> "$COMBINED"

            printf "%s\t%s\t%s\tREADY\n" \
                "$batch" "$RID" "$nq" >> "$SUMMARY"

            echo "Saved: $RESULT"

            success=1
            break
        fi

        echo "Could not parse status; retrying..."
        sleep "$POLL_INTERVAL"

    done

    [[ "$success" -eq 0 ]] && \
        echo "WARNING: $batch did not complete successfully."

    echo

done


echo "============================================================"
echo "BLASTx batches finished"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

echo
echo "Contigs with protein hits:"

if [[ -s "$COMBINED" ]]; then
    cut -d',' -f1 "$COMBINED" \
    | sed '/^[[:space:]]*$/d' \
    | sort -u
else
    echo "None"
fi
