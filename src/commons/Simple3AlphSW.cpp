#include "Simple3AlphSW.h"
#include <cstring>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>

Simple3AlphSW::Simple3AlphSW()
    : queryProfile(NULL),
      scoreMatrix(NULL),
      queryGapMatrix(NULL),
      targetGapMatrix(NULL),
      backtraceMatrix(NULL),
      queryLen(0),
      qSeqAA(NULL),
      qSeq3Di(NULL),
      qSeq12st(NULL),
      targetLen(0),
      tSeqAA(NULL),
      tSeq3Di(NULL),
      tSeq12st(NULL),
      maxQueryLen(0),
      maxTargetLen(0),
      maxScore(0),
      maxScoreRow(0),
      maxScoreCol(0),
      alignStartRow(0),
      alignStartCol(0),
      lastGapOpen(0),
      lastGapExtend(0)
{
    // Memory will be allocated in init_ssw based on actual query length
}

Simple3AlphSW::~Simple3AlphSW() {
    freeMemory();
}

void Simple3AlphSW::allocateMemory() {
    // Allocate query profile as flat array: [maxQueryLen * N_ALPHABETS * MAX_ALPH_SIZE]
    size_t queryProfileSize = maxQueryLen * N_ALPHABETS * MAX_ALPH_SIZE;
    queryProfile = new int16_t[queryProfileSize];
    std::memset(queryProfile, 0, queryProfileSize * sizeof(int16_t));

    // Allocate DP matrices as flat arrays: [(maxQueryLen+1) * (maxTargetLen+1)]
    size_t matrixSize = (maxQueryLen + 1) * (maxTargetLen + 1);

    scoreMatrix = new int16_t[matrixSize];
    std::memset(scoreMatrix, 0, matrixSize * sizeof(int16_t));

    queryGapMatrix = new int16_t[matrixSize];
    std::memset(queryGapMatrix, 0, matrixSize * sizeof(int16_t));

    targetGapMatrix = new int16_t[matrixSize];
    std::memset(targetGapMatrix, 0, matrixSize * sizeof(int16_t));

    backtraceMatrix = new uint8_t[matrixSize];
    std::memset(backtraceMatrix, BT_TERM, matrixSize * sizeof(uint8_t));
}

void Simple3AlphSW::freeMemory() {
    // Free flat arrays - simple single deletes
    if (queryProfile != NULL) {
        delete[] queryProfile;
        queryProfile = NULL;
    }

    if (scoreMatrix != NULL) {
        delete[] scoreMatrix;
        scoreMatrix = NULL;
    }

    if (queryGapMatrix != NULL) {
        delete[] queryGapMatrix;
        queryGapMatrix = NULL;
    }

    if (targetGapMatrix != NULL) {
        delete[] targetGapMatrix;
        targetGapMatrix = NULL;
    }

    if (backtraceMatrix != NULL) {
        delete[] backtraceMatrix;
        backtraceMatrix = NULL;
    }
}

void Simple3AlphSW::init_ssw(
    const unsigned char* q_aa_seq,
    const unsigned char* q_3di_seq,
    const unsigned char* q_12st_seq,
    int query_len,
    int8_t** submat_aa,
    int8_t** submat_3di,
    int8_t** submat_12st,
    size_t maxTargetLen_param)
{
    // Free existing memory if already allocated (BEFORE updating any variables)
    if (queryProfile != NULL) {
        freeMemory();
    }

    // Store query sequences (just pointers, no ownership)
    qSeqAA = q_aa_seq;
    qSeq3Di = q_3di_seq;
    qSeq12st = q_12st_seq;
    queryLen = query_len;

    // Set max target length (REQUIRED parameter)
    maxTargetLen = maxTargetLen_param;

    // Set max query length to actual query length
    maxQueryLen = queryLen;

    // Allocate memory based on actual query length
    allocateMemory();

    // Build query profile
    buildQueryProfile(submat_aa, submat_3di, submat_12st);
}

void Simple3AlphSW::buildQueryProfile(
    int8_t** submat_aa,
    int8_t** submat_3di,
    int8_t** submat_12st)
{
    // For each query position, precompute scores for all possible target letters
    // Access pattern: queryProfile[i * N_ALPHABETS * MAX_ALPH_SIZE + n * MAX_ALPH_SIZE + k]
    for (int i = 0; i < queryLen; i++) {
        int base_idx = i * N_ALPHABETS * MAX_ALPH_SIZE;

        // AA alphabet (21x21)
        unsigned char qaa = qSeqAA[i];
        int aa_offset = base_idx + 0 * MAX_ALPH_SIZE;
        for (int k = 0; k < SIZE_AA; k++) {
            queryProfile[aa_offset + k] = submat_aa[qaa][k];
        }

        // 3Di alphabet (21x21)
        unsigned char q3di = qSeq3Di[i];
        int di_offset = base_idx + 1 * MAX_ALPH_SIZE;
        for (int k = 0; k < SIZE_3DI; k++) {
            queryProfile[di_offset + k] = submat_3di[q3di][k];
        }

        // 12st alphabet (13x13)
        unsigned char q12st = qSeq12st[i];
        int st_offset = base_idx + 2 * MAX_ALPH_SIZE;
        for (int k = 0; k < SIZE_12ST; k++) {
            queryProfile[st_offset + k] = submat_12st[q12st][k];
        }
    }
}

int Simple3AlphSW::align(
    const unsigned char* t_aa_seq,
    const unsigned char* t_3di_seq,
    const unsigned char* t_12st_seq,
    int target_len,
    int gap_open,
    int gap_extend,
    bool compute_backtrace)
{
    // Store target sequences for backtrace
    tSeqAA = t_aa_seq;
    tSeq3Di = t_3di_seq;
    tSeq12st = t_12st_seq;
    targetLen = target_len;

    // Store gap penalties for backtrace
    lastGapOpen = gap_open;
    lastGapExtend = gap_extend;

    // Fill DP matrix (with optional backtrace)
    fillDPMatrix(t_aa_seq, t_3di_seq, t_12st_seq, target_len, gap_open, gap_extend, compute_backtrace);

    return maxScore;
}

void Simple3AlphSW::fillDPMatrix(
    const unsigned char* tSeqAA,
    const unsigned char* tSeq3Di,
    const unsigned char* tSeq12st,
    int targetLen,
    int gap_open,
    int gap_extend,
    bool compute_backtrace)
{
    // Matrix dimensions
    size_t rowSize = maxTargetLen + 1;

    // Initialize first row and column to 0 (local alignment)
    for (int i = 0; i <= queryLen; i++) {
        size_t idx = i * rowSize;  // i * rowSize + 0
        scoreMatrix[idx] = 0;
        queryGapMatrix[idx] = 0;
        targetGapMatrix[idx] = 0;
        if (compute_backtrace) {
            backtraceMatrix[idx] = BT_TERM;
        }
    }

    for (int j = 0; j <= targetLen; j++) {
        scoreMatrix[j] = 0;  // 0 * rowSize + j
        queryGapMatrix[j] = 0;
        targetGapMatrix[j] = 0;
        if (compute_backtrace) {
            backtraceMatrix[j] = BT_TERM;
        }
    }

    maxScore = 0;
    maxScoreRow = 0;
    maxScoreCol = 0;

    // Fill DP matrix using Smith-Waterman algorithm
    // Based on nalph_sw.py fill_dp_matrix function (lines 34-51)
    for (int i = 1; i <= queryLen; i++) {
        size_t row_idx = i * rowSize;
        size_t prev_row_idx = (i - 1) * rowSize;

        for (int j = 1; j <= targetLen; j++) {
            size_t idx = row_idx + j;
            size_t left_idx = row_idx + (j - 1);
            size_t up_idx = prev_row_idx + j;
            size_t diag_idx = prev_row_idx + (j - 1);

            // Update gap matrices
            // Query gap: gap in query (insertion in target, move left)
            int qgap_extend = queryGapMatrix[left_idx] - gap_extend;
            int qgap_open = scoreMatrix[left_idx] - gap_open;
            queryGapMatrix[idx] = std::max(qgap_extend, qgap_open);

            // Target gap: gap in target (deletion in target, move up)
            int tgap_extend = targetGapMatrix[up_idx] - gap_extend;
            int tgap_open = scoreMatrix[up_idx] - gap_open;
            targetGapMatrix[idx] = std::max(tgap_extend, tgap_open);

            // Sum scores from all 3 alphabets
            // Index i-1 because query profile is 0-indexed, but DP matrix is 1-indexed
            // Index j-1 because target sequences are 0-indexed
            int query_base = (i - 1) * N_ALPHABETS * MAX_ALPH_SIZE;
            int cellScore = 0;
            cellScore += queryProfile[query_base + 0 * MAX_ALPH_SIZE + tSeqAA[j - 1]];     // AA alphabet
            cellScore += queryProfile[query_base + 1 * MAX_ALPH_SIZE + tSeq3Di[j - 1]];    // 3Di alphabet
            cellScore += queryProfile[query_base + 2 * MAX_ALPH_SIZE + tSeq12st[j - 1]];   // 12st alphabet

            // Smith-Waterman recurrence (local alignment)
            int matchScore = scoreMatrix[diag_idx] + cellScore;

            // Take maximum of: 0 (local start), query gap, target gap, match
            // Order matches Python reference: max(0, q_gap, t_gap, match)
            // Use > to favor later options in ties: match > target gap > query gap > 0
            int maxVal = 0;
            uint8_t bt_state = BT_TERM;

            if (queryGapMatrix[idx] > maxVal) {
                maxVal = queryGapMatrix[idx];
                bt_state = BT_QGAP;
            }
            if (targetGapMatrix[idx] > maxVal) {
                maxVal = targetGapMatrix[idx];
                bt_state = BT_TGAP;
            }
            if (matchScore > maxVal) {
                maxVal = matchScore;
                bt_state = BT_MATCH;
            }

            scoreMatrix[idx] = maxVal;

            // Only store backtrace if requested
            if (compute_backtrace) {
                backtraceMatrix[idx] = bt_state;
            }

            // Track maximum score position
            // Use > to keep first occurrence (lexicographically smallest position)
            if (scoreMatrix[idx] > maxScore) {
                maxScore = scoreMatrix[idx];
                maxScoreRow = i;
                maxScoreCol = j;
            }
        }
    }
}

std::string Simple3AlphSW::getBacktrace() {
    std::string result;
    backtrace(result);
    return result;
}

void Simple3AlphSW::backtrace(std::string& result) {
    // Start from maximum score position and work backwards
    // Based on nalph_sw.py backtrace_alignment function (lines 224-266)
    // Uses state machine: S_MAT, QGAP_MAT, TGAP_MAT
    const int S_MAT = 0, QGAP_MAT = 1, TGAP_MAT = 2;

    int i = maxScoreRow;
    int j = maxScoreCol;
    int pos = S_MAT;  // Current state

    result.clear();

    size_t rowSize = maxTargetLen + 1;

    // Backtrace with state machine (matches Python exactly)
    // Python condition: while (s[i, j] > 0 or pos != S_MAT) and (i > 1 and j > 1)
    // But we need to process position 1 (sequence index 0), so use (i > 0 && j > 0)
    while ((scoreMatrix[i * rowSize + j] > 0 || pos != S_MAT) && (i > 0 && j > 0)) {
        size_t idx = i * rowSize + j;
        size_t left_idx = i * rowSize + (j - 1);
        size_t up_idx = (i - 1) * rowSize + j;
        size_t diag_idx = (i - 1) * rowSize + (j - 1);

        if (pos == S_MAT) {
            // Check which state to follow (Python checks: q_gap, t_gap, match)
            if (scoreMatrix[idx] == queryGapMatrix[idx]) {
                // Check if we're extending a gap
                if (queryGapMatrix[idx] == queryGapMatrix[left_idx] - lastGapExtend) {
                    pos = QGAP_MAT;
                }
                result += 'D';
                j--;
            } else if (scoreMatrix[idx] == targetGapMatrix[idx]) {
                // Check if we're extending a gap
                if (targetGapMatrix[idx] == targetGapMatrix[up_idx] - lastGapExtend) {
                    pos = TGAP_MAT;
                }
                result += 'I';
                i--;
            } else {
                // Match/mismatch: check if score came from diagonal
                int query_base = (i - 1) * N_ALPHABETS * MAX_ALPH_SIZE;
                int cellScore = 0;
                cellScore += queryProfile[query_base + 0 * MAX_ALPH_SIZE + tSeqAA[j - 1]];
                cellScore += queryProfile[query_base + 1 * MAX_ALPH_SIZE + tSeq3Di[j - 1]];
                cellScore += queryProfile[query_base + 2 * MAX_ALPH_SIZE + tSeq12st[j - 1]];

                if (scoreMatrix[idx] == scoreMatrix[diag_idx] + cellScore) {
                    result += 'M';
                    i--;
                    j--;
                } else {
                    // Shouldn't happen - score doesn't match any source
                    break;
                }
            }
        } else if (pos == QGAP_MAT) {
            // In query gap state: check if we're at gap open
            if (queryGapMatrix[idx] == scoreMatrix[left_idx] - lastGapOpen) {
                pos = S_MAT;
            }
            result += 'D';
            j--;
        } else if (pos == TGAP_MAT) {
            // In target gap state: check if we're at gap open
            if (targetGapMatrix[idx] == scoreMatrix[up_idx] - lastGapOpen) {
                pos = S_MAT;
            }
            result += 'I';
            i--;
        }
    }

    // Store alignment start positions
    // After loop, i and j are already at the position PAST the last processed cell
    // So they already represent 0-based sequence indices
    alignStartRow = i;
    alignStartCol = j;

    // Reverse the path (was built backwards)
    std::reverse(result.begin(), result.end());
}
