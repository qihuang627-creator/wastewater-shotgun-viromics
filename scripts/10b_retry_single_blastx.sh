#!/usr/bin/env bash

set -u -o pipefail

IDS="results_summary/blastx_retry_candidates.txt"
INDIR="results_summary/blastx_single"

OUTDIR="results/protein_blast"
LOGDIR="logs"

API="https://blast.ncbi.nlm.nih.gov/Blast.cgi"

COMBINED="${OUTDIR}/ERR14789190_blastx_single_retry.csv"
SUMMARY="${OUTDIR}/ERR14789190_blastx_single_retry_jobs.tsv"

POLL_INTERVAL=60
MAX_POLLS=180

mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$IDS" ]]; then
    echo "ERROR: Missing retry candidate list:"
    echo "$IDS"
    exit 1
fi

: > "$COMBINED"

printf "contig\tRID\tstatus\n" > "$SUMMARY"


while read -r CONTIG; do

    [[ -z "$CONTIG" ]] && continue

    QUERY="${INDIR}/${CONTIG}.fa"

    if [[ ! -s "$QUERY" ]]; then
        echo "ERROR: Missing FASTA: $QUERY"
        printf "%s\tNA\tMISSING_FASTA\n" \
            "$CONTIG" >> "$SUMMARY"
        continue
    fi

    SUBMIT="${OUTDIR}/${CONTIG}_blastx_submit.txt"
    STATUS_FILE="${OUTDIR}/.${CONTIG}_blastx_status.tmp"
    RESULT_TMP="${OUTDIR}/.${CONTIG}_blastx_result.tmp"
    RESULT="${OUTDIR}/${CONTIG}_blastx.csv"

    echo
    echo "============================================================"
    echo "Starting BLASTx: $CONTIG"
    echo "============================================================"

    # --------------------------------------------------------
    # Submit one contig
    # --------------------------------------------------------

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
        printf "%s\tNA\tSUBMISSION_FAILED\n" \
            "$CONTIG" >> "$SUMMARY"
        continue
    fi


    RID=$(
        awk -F'=' '
        /RID =/ {
            gsub(/[[:space:]]/, "", $2)
            print $2
            exit
        }
        ' "$SUBMIT"
    )

    RTOE=$(
        awk -F'=' '
        /RTOE =/ {
            gsub(/[[:space:]]/, "", $2)
            print $2
            exit
        }
        ' "$SUBMIT"
    )


    if [[ -z "$RID" ]]; then
        echo "ERROR: No RID returned."
        printf "%s\tNA\tNO_RID\n" \
            "$CONTIG" >> "$SUMMARY"
        continue
    fi


    echo "RID: $RID"

    [[ -z "$RTOE" ]] && RTOE=60

    if (( RTOE < 60 )); then
        WAIT=60
    else
        WAIT=$RTOE
    fi

    echo "Initial wait: ${WAIT}s"
    sleep "$WAIT"


    SUCCESS=0

    # --------------------------------------------------------
    # Poll
    # --------------------------------------------------------

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
            echo "Network polling failed; retrying..."
            sleep "$POLL_INTERVAL"
            continue
        fi


        if grep -qiE \
            'CPU usage limit was exceeded|SIGXCPU|An error has occurred' \
            "$STATUS_FILE"
        then
            echo "NCBI server-side error."
            printf "%s\t%s\tSERVER_ERROR\n" \
                "$CONTIG" "$RID" >> "$SUMMARY"
            break
        fi


        STATUS=$(
            grep -oE \
                'Status=(WAITING|READY|FAILED|UNKNOWN)' \
                "$STATUS_FILE" \
            | head -1 \
            | cut -d= -f2
        )


        case "$STATUS" in

            WAITING)

                echo "Status: WAITING"
                sleep "$POLL_INTERVAL"
                ;;


            FAILED|UNKNOWN)

                echo "Status: $STATUS"

                printf "%s\t%s\t%s\n" \
                    "$CONTIG" "$RID" "$STATUS" \
                    >> "$SUMMARY"

                break
                ;;


            READY)

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
                    echo "Download failed; retrying..."
                    sleep "$POLL_INTERVAL"
                    continue
                fi


                if grep -qiE \
                    'CPU usage limit was exceeded|SIGXCPU|An error has occurred' \
                    "$RESULT_TMP"
                then

                    echo "Result contains NCBI error."

                    printf "%s\t%s\tRESULT_ERROR\n" \
                        "$CONTIG" "$RID" >> "$SUMMARY"

                    break
                fi


                # Remove empty lines
                sed '/^[[:space:]]*$/d' \
                    "$RESULT_TMP" \
                    > "$RESULT"


                if [[ -s "$RESULT" ]]; then

                    cat "$RESULT" >> "$COMBINED"

                    echo "Protein hit(s) retrieved."

                    printf "%s\t%s\tREADY_WITH_HITS\n" \
                        "$CONTIG" "$RID" >> "$SUMMARY"

                else

                    echo "No significant protein hits."

                    printf "%s\t%s\tREADY_NO_HITS\n" \
                        "$CONTIG" "$RID" >> "$SUMMARY"

                fi


                rm -f "$RESULT_TMP"

                SUCCESS=1
                break
                ;;


            *)

                echo "Could not parse status; retrying..."
                sleep "$POLL_INTERVAL"
                ;;

        esac

    done


    if [[ "$SUCCESS" -eq 0 ]]; then
        echo "WARNING: $CONTIG did not complete successfully."
    fi

done < "$IDS"


echo
echo "============================================================"
echo "Single-contig BLASTx retries finished"
echo "============================================================"

column -t -s $'\t' "$SUMMARY"

echo
echo "Contigs with retrieved protein hits:"

if [[ -s "$COMBINED" ]]; then

    cut -d',' -f1 "$COMBINED" \
    | sed '/^[[:space:]]*$/d' \
    | sort -u

else
    echo "None"
fi
