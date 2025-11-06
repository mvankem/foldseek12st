#ifndef SIMPLE3ALPHSW_H
#define SIMPLE3ALPHSW_H

#include <string>
#include <algorithm>
#include <cstdint>

/**
 * Simple3AlphSW: Simple (non-SIMD) Smith-Waterman alignment for 3 alphabets
 *
 * Performs Smith-Waterman alignment using three structural/sequence alphabets:
 * - AA (amino acid): 21 letters
 * - 3Di (structural alphabet): 21 letters
 * - 12st (structural alphabet): 13 letters
 *
 * Based on nalph_sw.py implementation - straightforward DP with query profiles.
 * Scores from all alphabets are summed at each position.
 */
class Simple3AlphSW {
public:
    /**
     * Constructor: initialize pointers (memory allocated in init_ssw)
     */
    Simple3AlphSW();

    /**
     * Destructor: free all allocated memory
     */
    ~Simple3AlphSW();

    /**
     * Initialize with query sequences and substitution matrices
     * Allocates memory based on query length and builds query profile
     *
     * @param q_aa_seq Query AA numeric sequence (0-based indices)
     * @param q_3di_seq Query 3Di numeric sequence (0-based indices)
     * @param q_12st_seq Query 12st numeric sequence (0-based indices)
     * @param query_len Query length (same for all 3 alphabets)
     * @param submat_aa AA substitution matrix as int8_t** [21][21]
     * @param submat_3di 3Di substitution matrix as int8_t** [21][21]
     * @param submat_12st 12st substitution matrix as int8_t** [13][13]
     * @param maxTargetLen Maximum expected target length (REQUIRED)
     */
    void init_ssw(
        const unsigned char* q_aa_seq,
        const unsigned char* q_3di_seq,
        const unsigned char* q_12st_seq,
        int query_len,
        int8_t** submat_aa,
        int8_t** submat_3di,
        int8_t** submat_12st,
        size_t maxTargetLen
    );

    /**
     * Perform Smith-Waterman alignment against target sequences
     *
     * @param t_aa_seq Target AA numeric sequence (0-based indices)
     * @param t_3di_seq Target 3Di numeric sequence (0-based indices)
     * @param t_12st_seq Target 12st numeric sequence (0-based indices)
     * @param target_len Target length (same for all 3 alphabets)
     * @param gap_open Gap open penalty (positive value)
     * @param gap_extend Gap extend penalty (positive value)
     * @param compute_backtrace If true, store backtrace matrix for getBacktrace() (default: false)
     * @return Maximum alignment score
     */
    int align(
        const unsigned char* t_aa_seq,
        const unsigned char* t_3di_seq,
        const unsigned char* t_12st_seq,
        int target_len,
        int gap_open,
        int gap_extend,
        bool compute_backtrace = false
    );

    /**
     * Get backtrace alignment path from last align() call
     * NOTE: align() must be called with compute_backtrace=true first
     * @return String with alignment path: 'M' (match), 'I' (insertion), 'D' (deletion)
     */
    std::string getBacktrace();

    /**
     * Get maximum score from last alignment
     */
    int getMaxScore() const { return maxScore; }

    /**
     * Get position of maximum score (query coordinate)
     */
    int getMaxScoreRow() const { return maxScoreRow; }

    /**
     * Get position of maximum score (target coordinate)
     */
    int getMaxScoreCol() const { return maxScoreCol; }

    /**
     * Get alignment start position (query coordinate, 0-based)
     * NOTE: Only valid after getBacktrace() has been called
     */
    int getAlignStartRow() const { return alignStartRow; }

    /**
     * Get alignment start position (target coordinate, 0-based)
     * NOTE: Only valid after getBacktrace() has been called
     */
    int getAlignStartCol() const { return alignStartCol; }

private:
    // Alphabet sizes
    static const int SIZE_AA = 21;
    static const int SIZE_3DI = 21;
    static const int SIZE_12ST = 13;
    static const int MAX_ALPH_SIZE = 21;  // max(21, 21, 13)
    static const int N_ALPHABETS = 3;

    // Backtrace states
    static const uint8_t BT_MATCH = 0;   // Match/mismatch (diagonal)
    static const uint8_t BT_QGAP = 1;    // Query gap (insertion in target, move left)
    static const uint8_t BT_TGAP = 2;    // Target gap (deletion in target, move up)
    static const uint8_t BT_TERM = 255;  // Terminal (score = 0)

    // Query profile: flat array [query_len * N_ALPHABETS * MAX_ALPH_SIZE]
    // Access: queryProfile[i * N_ALPHABETS * MAX_ALPH_SIZE + n * MAX_ALPH_SIZE + k]
    // For position i, alphabet n, stores scores for all possible target chars k
    int16_t* queryProfile;

    // DP matrices: flat arrays [(maxQueryLen+1) * (maxTargetLen+1)]
    // Access: matrix[i * (maxTargetLen+1) + j]
    int16_t* scoreMatrix;       // Main scores (S)
    int16_t* queryGapMatrix;    // Query gap scores (Q_gap)
    int16_t* targetGapMatrix;   // Target gap scores (T_gap)

    // Backtrace matrix to reconstruct alignment
    uint8_t* backtraceMatrix;   // BT_MATCH, BT_QGAP, BT_TGAP, BT_TERM

    // Query data (pointers to mapped sequences)
    int queryLen;
    const unsigned char* qSeqAA;
    const unsigned char* qSeq3Di;
    const unsigned char* qSeq12st;

    // Target data (pointers from last align() call)
    int targetLen;
    const unsigned char* tSeqAA;
    const unsigned char* tSeq3Di;
    const unsigned char* tSeq12st;

    // Max dimensions
    size_t maxQueryLen;
    size_t maxTargetLen;

    // Alignment result (from last align() call)
    int maxScore;
    int maxScoreRow;
    int maxScoreCol;

    // Alignment start positions (from last backtrace() call, 0-based)
    int alignStartRow;
    int alignStartCol;

    // Helper methods

    /**
     * Allocate all memory for DP matrices and query profile
     */
    void allocateMemory();

    /**
     * Free all allocated memory
     */
    void freeMemory();

    /**
     * Build query profile from query sequences and substitution matrices
     * For each query position and alphabet, precompute scores for all possible target letters
     */
    void buildQueryProfile(
        int8_t** submat_aa,
        int8_t** submat_3di,
        int8_t** submat_12st
    );

    /**
     * Fill DP matrices using Smith-Waterman algorithm with 3 alphabets
     * Scores from all alphabets are summed at each position
     * @param compute_backtrace If true, store backtrace matrix for path reconstruction
     */
    void fillDPMatrix(
        const unsigned char* tSeqAA,
        const unsigned char* tSeq3Di,
        const unsigned char* tSeq12st,
        int targetLen,
        int gap_open,
        int gap_extend,
        bool compute_backtrace
    );

    /**
     * Backtrace from maximum score position to build alignment path
     */
    void backtrace(std::string& result);
};

#endif // SIMPLE3ALPHSW_H
