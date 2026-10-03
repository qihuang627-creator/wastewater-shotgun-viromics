#!/usr/bin/env bash

set -u -o pipefail

# ============================================================
# Configuration
# ============================================================

IDS="results_summary/sensitive_nt_candidates.txt"
SOURCE_FA="results_summary/high_priority_panel_external_contigs.fa"

BATCHDIR="results_summary/sensitive_core_nt_batches"
OUTDIR="results/broad_blast"
LOGDIR="logs"

API="https://blast.ncbi.nlm.nih.gov/Blast.cgi"

COMBINED="${OUTDIR}/ERR14789190_sensitive_core_nt_blastn.csv"
SUMMARY="${OUTDIR}/ERR14789190_sensitive_core_nt_jobs.tsv"

BATCH_SIZE=5
POLL_INTERVAL=60
MAX_POLLS=180

mkdir -p "$BATCHDIR" "$OUTDIR" "$LOGDIR"


# ============================================================
# Check inputs
# ============================================================

if [[ ! -s "$IDS" ]]; then
    echo "ERROR: Missing candidate list:"
    echo "$IDS"
    exit 1
fi

if [[ ! -s "$SOURCE_FA" ]]; then
    echo "ERROR: Missing source FASTA:"
    echo "$SOURCE_FA"
    exit 1
fi

N_TOTAL=$(wc -l < "$IDS")

echo "Total candidate contigs: $N_TOTAL"
echo "Batch size: $BATCH_SIZE"
echo "Database: core_nt"
echo "Program: standard blastn"
echo "Word size: 11"
echo


# ============================================================
# Prepare batches
# ============================================================

rm -f "${BATCHDIR}"/batch_*.ids
rm -f "${BATCHDIR}"/batch_*.fa
rm -f "${BATCHDIR}"/batch_*.csv

split \
    -l "$BATCH_SIZE" \
    -d \
    -a 2 \
    "$IDS" \
    "${BATCHDIR}/batch_"

for ids_file in "${BATCHDIR}"/batch_*; do

    [[ "$ids_file" == *.fa ]] && continue
    [[ "$ids_file" == *.csv ]] && continue

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


echo "Prepared batches:"

for fa in "${BATCHDIR}"/batch_*.fa; do
    echo "$(basename "$fa"): $(grep -c '^>' "$fa") sequences"
done

echo


# ============================================================
# Initialize output files
# ============================================================

: > "$COMBINED"

printf "batch\tRID\tn_queries\tstatus\n" > "$SUMMARY"


# ============================================================
# Process each batch sequentially
# ============================================================

for QUERY in "${BATCHDIR}"/batch_*.fa; do

    batch=$(basename "$QUERY" .fa)
    nq=$(grep -c "^>" "$QUERY")

    SUBMIT_RESPONSE="${OUTDIR}/${batch}_submit.txt"
    STATUS_FILE="${OUTDIR}/.${batch}_status.tmp"
    RESULT_TMP="${OUTDIR}/.${batch}_result.tmp"
    RESULT="${BATCHDIR}/${batch}.csv"

    echo "============================================================"
    echo "Starting $batch"
    echo "Queries: $nq"
    echo "============================================================"


    # --------------------------------------------------------
    # Submit
    # --------------------------------------------------------

    if ! curl -sS \
        --retry 3 \
        --retry-delay 10 \
        --connect-timeout 60 \
        --max-time 180 \
        --data-urlencode "CMD=Put" \
        --data-urlencode "PROGRAM=blastn" \
        --data-urlencode "DATABASE=core_nt" \
        --data-urlencode "WORD_SIZE=11" \
        --data-urlencode "EXPECT=1e-5" \
        --data-urlencode "HITLIST_SIZE=20" \
        --data-urlencode "FILTER=mL" \
        --data-urlencode "TOOL=wwshotgun_viromics" \
        --data-urlencode "QUERY@${QUERY}" \
        "$API" \
        > "$SUBMIT_RESPONSE"
    then
        echo "ERROR: Submission failed for $batch"
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
        echo "ERROR: No RID returned for $batch"
        printf "%s\tNA\t%s\tNO_RID\n" \
            "$batch" "$nq" >> "$SUMMARY"
        continue
    fi


    echo "RID: $RID"

    if [[ -n "$RTOE" ]]; then
        echo "NCBI estimated runtime: ${RTOE} seconds"
    else
        RTOE=60
    fi


    # NCBI recommends avoiding aggressive polling.
    if (( RTOE < 60 )); then
        WAIT=60
    else
        WAIT=$RTOE
    fi

    echo "Initial wait: ${WAIT} seconds"
    sleep "$WAIT"


    # --------------------------------------------------------
    # Poll SearchInfo
    # --------------------------------------------------------

    success=0

    for ((poll=1; poll<=MAX_POLLS; poll++)); do

        echo "$(date) Checking $batch RID $RID ..."

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
            echo "Network polling failed; retrying in ${POLL_INTERVAL}s..."
            sleep "$POLL_INTERVAL"
            continue
        fi


        # Explicit server-side error detection
        if grep -qiE \
            'CPU usage limit was exceeded|An error has occurred|SIGXCPU' \
            "$STATUS_FILE"
        then
            echo "ERROR: NCBI server-side failure for $batch"
            printf "%s\t%s\t%s\tSERVER_ERROR\n" \
                "$batch" "$RID" "$nq" >> "$SUMMARY"
            break
        fi


        STATUS=$(
            grep -oE \
                'Status=(WAITING|READY|FAILED|UNKNOWN)' \
                "$STATUS_FILE" \
            | head -n 1 \
            | cut -d'=' -f2
        )


        case "$STATUS" in

            WAITING)
                echo "Status: WAITING"
                sleep "$POLL_INTERVAL"
                ;;


            FAILED)
                echo "ERROR: NCBI reports FAILED"
                printf "%s\t%s\t%s\tFAILED\n" \
                    "$batch" "$RID" "$nq" >> "$SUMMARY"
                break
                ;;


            UNKNOWN)
                echo "ERROR: RID unknown or expired"
                printf "%s\t%s\t%s\tUNKNOWN\n" \
                    "$batch" "$RID" "$nq" >> "$SUMMARY"
                break
                ;;


            READY)

                echo "Status: READY"

                HITS=$(
                    grep -oE \
                        'ThereAreHits=(yes|no)' \
                        "$STATUS_FILE" \
                    | head -n 1 \
                    | cut -d'=' -f2
                )

                if [[ "$HITS" == "no" ]]; then

                    echo "No significant hits for $batch"

                    : > "$RESULT"

                    printf "%s\t%s\t%s\tREADY_NO_HITS\n" \
                        "$batch" "$RID" "$nq" >> "$SUMMARY"

                    success=1
                    break
                fi


                # --------------------------------------------
                # Retrieve CSV results
                # --------------------------------------------

                echo "Downloading result..."

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
                    echo "Result download failed; retrying..."
                    sleep "$POLL_INTERVAL"
                    continue
                fi


                # Reject error pages masquerading as results
                if grep -qiE \
                    'CPU usage limit was exceeded|An error has occurred|SIGXCPU|Status=FAILED|Status=UNKNOWN' \
                    "$RESULT_TMP"
                then
                    echo "ERROR: Result contains NCBI error message"

                    printf "%s\t%s\t%s\tRESULT_ERROR\n" \
                        "$batch" "$RID" "$nq" >> "$SUMMARY"

                    break
                fi


                mv "$RESULT_TMP" "$RESULT"

                cat "$RESULT" >> "$COMBINED"

                printf "%s\t%s\t%s\tREADY_WITH_HITS\n" \
                    "$batch" "$RID" "$nq" >> "$SUMMARY"

                echo "Saved: $RESULT"

                success=1
                break
                ;;


            *)
                echo "Could not parse status; retrying in ${POLL_INTERVAL}s..."
                sleep "$POLL_INTERVAL"
                ;;

        esac

    done


    if [[ "$success" -eq 0 ]]; then
        echo "WARNING: $batch did not complete successfully."
    fi

    echo

done


# ============================================================
# Final summary
# ============================================================

echo "============================================================"
echo "All batches processed"
echo "============================================================"

echo
echo "Job summary:"
column -t -s $'\t' "$SUMMARY"

echo
echo "Combined result:"
echo "$COMBINED"

if [[ -s "$COMBINED" ]]; then

    echo
    echo "Number of contigs with at least one nt hit:"

    cut -d',' -f1 "$COMBINED" \
    | sort -u \
    | wc -l

else

    echo
    echo "Combined BLAST result is empty."

fi
