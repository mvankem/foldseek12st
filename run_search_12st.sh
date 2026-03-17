#!/bin/sh -ex
QUERY="${DATADIR}/scop"
QUERYDB="${RESULTS}/query"
"${FOLDSEEK}" createdb "${QUERY}" "${QUERYDB}"

TARGET="${DATADIR}/scop"
SCOPANOTATION="${DATADIR}/scop_lookup_bench.tsv"
TARGETDB="${RESULTS}/target"
"${FOLDSEEK}" createdb "${TARGET}" "${TARGETDB}" 

"${FOLDSEEK}" createindex "$TARGETDB"  "$RESULTS/tmp" 


"${FOLDSEEK}" prefilter "${QUERYDB}_ss" "${TARGETDB}_ss" "$RESULTS/tmp/prefDB" -s 9 --max-seqs 100000 --comp-bias-corr 1 --comp-bias-corr-scale 0.15
"${FOLDSEEK}" rescorediagonal12st "$QUERYDB" "$TARGETDB" "$RESULTS/tmp/prefDB" "$RESULTS/tmp/prefDB_rescored"
"${FOLDSEEK}" structurealign "$QUERYDB" "$TARGETDB" "$RESULTS/tmp/prefDB_rescored" "$RESULTS/results_aln" --submat-12st-scale 2.1 --alignment-type 4 -e 10000 --max-accept 100
"${FOLDSEEK}" convertalis "$QUERYDB" "$TARGETDB" "$RESULTS/results_aln" "$RESULTS/results_aln.m8"

"${EVALUATE}" "$SCOPANOTATION" "$RESULTS/results_aln.m8" > "${RESULTS}/evaluation.log"

ACTUAL=$(awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' "${RESULTS}/evaluation.log")
TARGET="0.986667 0.796061 0.459869"
awk -v actual="$ACTUAL" -v target="$TARGET" \
    'BEGIN { print (actual >= target) ? "GOOD" : "BAD"; print "Expected: ", target; print "Actual: ", actual; }' \
    > "${RESULTS}.report"

# Expected:  0.986667 0.796061 0.459869
# Actual:  0.986667 0.833838 0.549373
