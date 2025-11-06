#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include "FileUtil.h"
#include "Util.h"
#include "Simple3AlphSW.h"
#include <iomanip>

const char* binary_name = "test_simple3alphsw";

// Alphabet encodings
const char* ALPHABET_AA = "ACDEFGHIKLMNPQRSTVWYX";  // 21 letters (AA and 3Di)
const char* ALPHABET_12ST = "ABCDEFGHIJKLX";         // 13 letters (12st)

// Helper function: encode a string sequence to numeric array using alphabet
// Returns allocated array (caller must delete[])
// Exits with error if unknown character is found
unsigned char* encodeSequence(const std::string& seq, const char* alphabet) {
    unsigned char* encoded = new unsigned char[seq.length()];

    for (size_t i = 0; i < seq.length(); i++) {
        char c = seq[i];
        // Find position in alphabet
        const char* pos = strchr(alphabet, c);
        if (pos != NULL) {
            encoded[i] = (unsigned char)(pos - alphabet);
        } else {
            // Unknown character - error out
            std::cerr << "Error: Unknown character '" << c << "' at position " << i
                      << " in sequence. Not found in alphabet: " << alphabet << std::endl;
            delete[] encoded;
            exit(1);
        }
    }

    return encoded;
}

// Helper function: convert backtrace string (M, I, D) to CIGAR format
std::string backtraceToCigar(const std::string& backtrace) {
    if (backtrace.empty()) {
        return "";
    }

    std::string cigar;
    char currentOp = backtrace[0];
    int count = 1;

    for (size_t i = 1; i < backtrace.length(); i++) {
        if (backtrace[i] == currentOp) {
            count++;
        } else {
            cigar += std::to_string(count) + currentOp;
            currentOp = backtrace[i];
            count = 1;
        }
    }
    // Add the last run
    cigar += std::to_string(count) + currentOp;

    return cigar;
}

// Helper function: generate aligned sequences from backtrace
void getAlignedSequences(const std::string& backtrace,
                        const std::string& q_seq,
                        const std::string& t_seq,
                        int q_start, int t_start,
                        std::string& q_aligned,
                        std::string& t_aligned) {
    q_aligned.clear();
    t_aligned.clear();

    int q_pos = q_start;
    int t_pos = t_start;

    for (char op : backtrace) {
        if (op == 'M') {
            q_aligned += q_seq[q_pos];
            t_aligned += t_seq[t_pos];
            q_pos++;
            t_pos++;
        } else if (op == 'I') {
            q_aligned += q_seq[q_pos];
            t_aligned += '-';
            q_pos++;
        } else if (op == 'D') {
            q_aligned += '-';
            t_aligned += t_seq[t_pos];
            t_pos++;
        }
    }
}

// Substitution matrices as float arrays (original values from comments)
// 3Di matrix: 21x21 (ACDEFGHIKLMNPQRSTVWYX) - added X row/column with zeros
const float SUBMAT_3DI_FLOAT[21][21] = {
    {2.811, -1.499, 0.495, 1.229, 1.408, -1.068, -0.803, -3.538, -1.683, -1.456, -5.126, -2.371, -0.596, 0.746, -1.984, -3.477, -2.323, -3.05, -0.069, -1.214, 0.0},
    {-1.499, 3.164, -0.953, -4.01, -2.357, -2.125, -1.904, -5.858, -6.345, 0.323, -6.906, 0.241, -0.07, 0.507, -0.258, 0.064, -3.794, 0.395, -3.686, -4.282, 0.0},
    {0.495, -0.953, 2.019, -1.657, -0.141, 0.269, 0.445, -1.668, -2.133, -2.177, -2.558, -1.105, 0.592, -0.65, -0.701, -2.069, -0.782, -1.672, -0.766, -1.223, 0.0},
    {1.229, -4.01, -1.657, 4.384, -0.995, -3.563, -1.934, -6.05, -5.028, -3.721, -8.559, -4.232, -3.067, -1.654, -3.853, -5.118, -4.753, -6.299, -3.097, -1.281, 0.0},
    {1.408, -2.357, -0.141, -0.995, 3.42, -1.323, -1.397, -2.5, 0.316, -1.395, -4.251, -2.651, -0.987, 0.957, -2.337, -3.823, -1.674, -3.374, 2.14, -2.235, 0.0},
    {-1.068, -2.125, 0.269, -3.563, -1.323, 3.215, 1.439, 0.155, -3.531, -3.539, -0.489, -0.842, -0.851, -1.794, 1.292, -1.436, 1.961, -2.986, -1.95, -0.993, 0.0},
    {-0.803, -1.904, 0.445, -1.934, -1.397, 1.439, 3.085, -1.961, -3.748, -3.062, -2.776, 0.157, -0.512, -1.527, 0.49, -1.649, -0.411, -2.581, -2.409, 1.354, 0.0},
    {-3.538, -5.858, -1.668, -6.05, -2.5, 0.155, -1.961, 4.149, -2.498, -5.65, 3.5, -3.415, -3.029, -2.995, -1.295, -4.439, 2.908, -5.954, -2.359, -3.802, 0.0},
    {-1.683, -6.345, -2.133, -5.028, 0.316, -3.531, -3.748, -2.498, 4.533, -5.443, -4.072, -5.853, -3.635, -2.577, -4.361, -7.156, -2.525, -7.332, 2.661, -3.777, 0.0},
    {-1.456, 0.323, -2.177, -3.721, -1.395, -3.539, -3.062, -5.65, -5.443, 2.911, -7.99, -1.313, -0.92, 1.017, -1.859, -1.781, -4.58, -0.05, -3.946, -4.496, 0.0},
    {-5.126, -6.906, -2.558, -8.559, -4.251, -0.489, -2.776, 3.5, -4.072, -7.99, 5.22, -4.701, -4.573, -5.149, -2.37, -5.17, 1.426, -7.763, -3.237, -4.54, 0.0},
    {-2.371, 0.241, -1.105, -4.232, -2.651, -0.842, 0.157, -3.415, -5.853, -1.313, -4.701, 3.251, 0.058, -1.036, 1.108, 1.385, -2.163, 0.036, -4.042, -2.268, 0.0},
    {-0.596, -0.07, 0.592, -3.067, -0.987, -0.851, -0.512, -3.029, -3.635, -0.92, -4.573, 0.058, 2.287, -0.114, -0.042, -0.937, -2.067, 0.143, -2.216, -2.938, 0.0},
    {0.746, 0.507, -0.65, -1.654, 0.957, -1.794, -1.527, -2.995, -2.577, 1.017, -5.149, -1.036, -0.114, 2.462, -1.21, -1.762, -2.532, -0.624, -0.95, -2.444, 0.0},
    {-1.984, -0.258, -0.701, -3.853, -2.337, 1.292, 0.49, -1.295, -4.361, -1.859, -2.37, 1.108, -0.042, -1.21, 2.832, 1.13, 0.113, -0.563, -3.08, -1.568, 0.0},
    {-3.477, 0.064, -2.069, -5.118, -3.823, -1.436, -1.649, -4.439, -7.156, -1.781, -5.17, 1.385, -0.937, -1.762, 1.13, 2.951, -2.836, 0.018, -5.458, -4.29, 0.0},
    {-2.323, -3.794, -0.782, -4.753, -1.674, 1.961, -0.411, 2.908, -2.525, -4.58, 1.426, -2.163, -2.067, -2.532, 0.113, -2.836, 3.985, -4.598, -2.288, -2.388, 0.0},
    {-3.05, 0.395, -1.672, -6.299, -3.374, -2.986, -2.581, -5.954, -7.332, -0.05, -7.763, 0.036, 0.143, -0.624, -0.563, 0.018, -4.598, 1.731, -5.092, -5.741, 0.0},
    {-0.069, -3.686, -0.766, -3.097, 2.14, -1.95, -2.409, -2.359, 2.661, -3.946, -3.237, -4.042, -2.216, -0.95, -3.08, -5.458, -2.288, -5.092, 3.979, -3.237, 0.0},
    {-1.214, -4.282, -1.223, -1.281, -2.235, -0.993, 1.354, -3.802, -3.777, -4.496, -4.54, -2.268, -2.938, -2.444, -1.568, -4.29, -2.388, -5.741, -3.237, 4.285, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}  // X row
};

// AA matrix: 21x21 (ACDEFGHIKLMNPQRSTVWYX)
const float SUBMAT_AA_FLOAT[21][21] = {
    {3.9291, -0.4085, -1.7534, -0.8639, -2.2101, 0.1596, -1.6251, -1.3218, -0.7340, -1.4646, -0.9353, -1.5307, -0.8143, -0.8040, -1.4135, 1.1158, -0.0454, -0.1894, -2.5269, -1.7640, -1.0000},
    {-0.4085, 8.5821, -3.4600, -3.6125, -2.3755, -2.5004, -2.9878, -1.2277, -3.0363, -1.2775, -1.4198, -2.6598, -2.7952, -2.9019, -3.3892, -0.8750, -0.8667, -0.8077, -2.3041, -2.4071, -1.0000},
    {-1.7534, -3.4600, 5.7742, 1.5103, -3.4839, -1.3135, -1.1189, -3.1212, -0.7018, -3.6057, -3.0585, 1.2717, -1.4801, -0.3134, -1.6058, -0.2610, -1.0507, -3.1426, -4.2143, -3.0650, -1.0000},
    {-0.8639, -3.6125, 1.5103, 4.9028, -3.1924, -2.1102, -0.1177, -3.1944, 0.7753, -2.8465, -1.9980, -0.2680, -1.1162, 1.8546, -0.1154, -0.1469, -0.8633, -2.4423, -2.8354, -2.0205, -1.0000},
    {-2.2101, -2.3755, -3.4839, -3.1924, 6.0461, -3.1074, -1.2342, -0.1609, -3.0787, 0.4148, 0.0126, -2.9940, -3.5973, -3.1644, -2.7863, -2.3690, -2.1076, -0.8490, 0.9176, 2.9391, -1.0000},
    {0.1596, -2.5004, -1.3135, -2.1102, -3.1074, 5.5633, -2.0409, -3.7249, -1.5280, -3.6270, -2.6766, -0.4228, -2.1335, -1.7852, -2.3041, -0.2925, -1.5754, -3.1387, -2.4915, -3.0398, -1.0000},
    {-1.6251, -2.9878, -1.1189, -0.1177, -1.2342, -2.0409, 7.5111, -3.2316, -0.7210, -2.7867, -1.5513, 0.5785, -2.1609, 0.4480, -0.2499, -0.8816, -1.6859, -3.1175, -2.3422, 1.6926, -1.0000},
    {-1.3218, -1.2277, -3.1212, -3.1944, -0.1609, -3.7249, -3.2316, 3.9985, -2.6701, 1.5216, 1.1268, -3.2170, -2.7567, -2.7696, -2.9902, -2.3482, -0.7176, 2.5470, -2.5805, -1.3314, -1.0000},
    {-0.7340, -3.0363, -0.7018, 0.7753, -3.0787, -1.5280, -0.7210, -2.6701, 4.5046, -2.4468, -1.3547, -0.1790, -1.0136, 1.2726, 2.1087, -0.2034, -0.6696, -2.2624, -2.9564, -1.8200, -1.0000},
    {-1.4646, -1.2775, -3.6057, -2.8465, 0.4148, -3.6270, -2.7867, 1.5216, -2.4468, 3.8494, 1.9918, -3.3789, -2.8601, -2.1339, -2.1546, -2.4426, -1.1975, 0.7884, -1.6319, -1.0621, -1.0000},
    {-0.9353, -1.4198, -3.0585, -1.9980, 0.0126, -2.6766, -1.5513, 1.1268, -1.3547, 1.9918, 5.3926, -2.1509, -2.4764, -0.4210, -1.3671, -1.4809, -0.6663, 0.6872, -1.4248, -0.9949, -1.0000},
    {-1.5307, -2.6598, 1.2717, -0.2680, -2.9940, -0.4228, 0.5785, -3.2170, -0.1790, -3.3789, -2.1509, 5.6532, -2.0004, 0.0017, -0.4398, 0.6009, -0.0461, -2.8763, -3.6959, -2.0818, -1.0000},
    {-0.8143, -2.7952, -1.4801, -1.1162, -3.5973, -2.1335, -2.1609, -2.7567, -1.0136, -2.8601, -2.4764, -2.0004, 7.3646, -1.2819, -2.1086, -0.8090, -1.0753, -2.3487, -3.6542, -2.9198, -1.0000},
    {-0.8040, -2.9019, -0.3134, 1.8546, -3.1644, -1.7852, 0.4480, -2.7696, 1.2726, -2.1339, -0.4210, 0.0017, -1.2819, 5.2851, 0.9828, -0.1011, -0.6753, -2.1984, -1.9465, -1.4211, -1.0000},
    {-1.4135, -3.3892, -1.6058, -0.1154, -2.7863, -2.3041, -0.2499, -2.9902, 2.1087, -2.1546, -1.3671, -0.4398, -2.1086, 0.9828, 5.4735, -0.7648, -1.1223, -2.5026, -2.6794, -1.6939, -1.0000},
    {1.1158, -0.8750, -0.2610, -0.1469, -2.3690, -0.2925, -0.8816, -2.3482, -0.2034, -2.4426, -1.4809, 0.6009, -0.8090, -0.1011, -0.7648, 3.8844, 1.3811, -1.6462, -2.7519, -1.6858, -1.0000},
    {-0.0454, -0.8667, -1.0507, -0.8633, -2.1076, -1.5754, -1.6859, -0.7176, -0.6696, -1.1975, -0.6663, -0.0461, -1.0753, -0.6753, -1.1223, 1.3811, 4.5453, -0.0555, -2.4289, -1.6060, -1.0000},
    {-0.1894, -0.8077, -3.1426, -2.4423, -0.8490, -3.1387, -3.1175, 2.5470, -2.2624, 0.7884, 0.6872, -2.8763, -2.3487, -2.1984, -2.5026, -1.6462, -0.0555, 3.7689, -2.8343, -1.2075, -1.0000},
    {-2.5269, -2.3041, -4.2143, -2.8354, 0.9176, -2.4915, -2.3422, -2.5805, -2.9564, -1.6319, -1.4248, -3.6959, -3.6542, -1.9465, -2.6794, -2.7519, -2.4289, -2.8343, 10.5040, 2.1542, -1.0000},
    {-1.7640, -2.4071, -3.0650, -2.0205, 2.9391, -3.0398, 1.6926, -1.3314, -1.8200, -1.0621, -0.9949, -2.0818, -2.9198, -1.4211, -1.6939, -1.6858, -1.6060, -1.2075, 2.1542, 6.5950, -1.0000},
    {-1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}
};

// 12st matrix: 13x13 (ABCDEFGHIJKLX)
const float SUBMAT_12ST_FLOAT[13][13] = {
    {2.431, -0.264, -2.745, -0.281, -2.328, -2.291, 0.02, 0.47, -2.831, -3.038, -0.462, -0.509, 0.358},
    {-0.264, 2.413, -2.581, -1.033, 0.289, -2.353, -0.856, -0.108, -1.661, -2.613, -0.194, -0.896, -0.995},
    {-2.745, -2.581, 2.854, 1.199, -1.569, -0.289, -1.433, -2.635, -2.051, -2.914, -3.225, -3.323, -0.762},
    {-0.281, -1.033, 1.199, 3.193, -2.111, -1.752, -0.151, -0.191, -3.24, -3.594, -1.911, -1.778, 0.664},
    {-2.328, 0.289, -1.569, -2.111, 2.172, -1.001, -1.827, -1.929, 0.057, -1.653, -1.331, -1.952, -1.495},
    {-2.291, -2.353, -0.289, -1.752, -1.001, 1.935, 0.202, -2.28, -1.835, -1.885, -2.793, -2.311, -1.689},
    {0.02, -0.856, -1.433, -0.151, -1.827, 0.202, 2.239, 0.002, -2.593, -2.68, -1.425, -1.284, -0.453},
    {0.47, -0.108, -2.635, -0.191, -1.929, -2.28, 0.002, 2.437, -2.416, -3.041, -0.013, -0.495, -0.113},
    {-2.831, -1.661, -2.051, -3.24, 0.057, -1.835, -2.593, -2.416, 2.647, 0.097, 0.787, -1.1, -0.75},
    {-3.038, -2.613, -2.914, -3.594, -1.653, -1.885, -2.68, -3.041, 0.097, 3.84, -1.298, 2.202, -0.117},
    {-0.462, -0.194, -3.225, -1.911, -1.331, -2.793, -1.425, -0.013, 0.787, -1.298, 2.954, 0.347, 0.25},
    {-0.509, -0.896, -3.323, -1.778, -1.952, -2.311, -1.284, -0.495, -1.1, 2.202, 0.347, 4.16, 0.72},
    {0.358, -0.995, -0.762, 0.664, -1.495, -1.689, -0.453, -0.113, -0.75, -0.117, 0.25, 0.72, 5.802}
};

// Helper function to apply scaling and discretize float matrix to int8_t
// Returns int8_t** pointer array (caller must delete[] both the pointers and data)
template<int ROWS, int COLS>
int8_t** scaleAndDiscretizeSubmat(const float floatData[ROWS][COLS], float scaleFactor) {
    // Allocate int8_t matrix
    int8_t** discretized = new int8_t*[ROWS];
    for (int i = 0; i < ROWS; i++) {
        discretized[i] = new int8_t[COLS];
        for (int j = 0; j < COLS; j++) {
            // Apply scaling and round to int8_t
            float scaled = floatData[i][j] * scaleFactor;
            discretized[i][j] = (int8_t)round(scaled);
        }
    }
    return discretized;
}

// Helper function to free int8_t** matrix
template<int ROWS>
void freeSubmat(int8_t** submat) {
    for (int i = 0; i < ROWS; i++) {
        delete[] submat[i];
    }
    delete[] submat;
}

int main(int argc, const char** argv) {
    // Open test data file
    const char* testDataFile = "sw_test_data.tsv";
    FILE* fp = FileUtil::openFileOrDie(testDataFile, "r", true);

    // Submat scaling factors (defined once for all tests)
    const float SUBMAT_AA_SCALE = 1.4f / 2.0f;   // = 0.7
    const float SUBMAT_3DI_SCALE = 2.1f;
    const float SUBMAT_12ST_SCALE = 2.1f;

    // Apply scaling and discretize substitution matrices (once for all tests)
    int8_t** submat_aa = scaleAndDiscretizeSubmat<21, 21>(SUBMAT_AA_FLOAT, SUBMAT_AA_SCALE);
    int8_t** submat_3di = scaleAndDiscretizeSubmat<21, 21>(SUBMAT_3DI_FLOAT, SUBMAT_3DI_SCALE);
    int8_t** submat_12st = scaleAndDiscretizeSubmat<13, 13>(SUBMAT_12ST_FLOAT, SUBMAT_12ST_SCALE);

    // // Set submat_aa to zero
    // for (int i = 0; i < 21; i++) {
    //     for (int j = 0; j < 21; j++) {
    //         submat_aa[i][j] = 0;
    //     }
    // }

    // // Set submat_12st to zero
    // for (int i = 0; i < 13; i++) {
    //     for (int j = 0; j < 13; j++) {
    //         submat_12st[i][j] = 0;
    //     }
    // }

    // // Print submat_3di
    // std::cout << "Substitution Matrix 3Di (scaled and discretized):" << std::endl;
    // for (int i = 0; i < 21; i++) {
    //     for (int j = 0; j < 21; j++) {
    //         // Pad ints for better readability
    //         std::cout << std::setw(4) << (int)submat_3di[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Track test results
    int totalTests = 0;
    int passedTests = 0;
    int failedTests = 0;

    // Create Simple3AlphSW instance (reused for all tests)
    Simple3AlphSW aligner;

    // Read all lines
    char buffer[100000];  // Large buffer for long sequences
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        totalTests++;

        //if (totalTests != 1) continue;

        // Parse tab-separated fields
        // Format: qname, tname, score, cigar, q_start, t_start, q_3di_seq, t_3di_seq, q_aa_seq, t_aa_seq, q_12st_seq, t_12st_seq
        const char* fields[15];
        size_t nFields = Util::getFieldsOfLine(buffer, fields, 15);

        if (nFields < 12) {
            std::cerr << "Error: Line " << totalTests << " - Expected 12 fields, got " << nFields << std::endl;
            failedTests++;
            continue;
        }

        // Extract field lengths (up to next tab or end of string)
        size_t qnameLen = Util::skipNonTab(fields[0]);
        size_t tnameLen = Util::skipNonTab(fields[1]);
        size_t scoreStrLen = Util::skipNonTab(fields[2]);
        size_t cigarLen = Util::skipNonTab(fields[3]);
        size_t qStartStrLen = Util::skipNonTab(fields[4]);
        size_t tStartStrLen = Util::skipNonTab(fields[5]);
        size_t q3diSeqLen = Util::skipNonTab(fields[6]);   // 3Di is first
        size_t t3diSeqLen = Util::skipNonTab(fields[7]);
        size_t qAaSeqLen = Util::skipNonTab(fields[8]);    // AA is second
        size_t tAaSeqLen = Util::skipNonTab(fields[9]);
        size_t q12stSeqLen = Util::skipNonTab(fields[10]); // 12st is third
        size_t t12stSeqLen = Util::skipNonTab(fields[11]);

        // Create strings from fields
        std::string qname(fields[0], qnameLen);
        std::string tname(fields[1], tnameLen);
        std::string scoreStr(fields[2], scoreStrLen);
        std::string expectedCigar(fields[3], cigarLen);
        std::string qStartStr(fields[4], qStartStrLen);
        std::string tStartStr(fields[5], tStartStrLen);
        std::string q_3di_seq(fields[6], q3diSeqLen);
        std::string t_3di_seq(fields[7], t3diSeqLen);
        std::string q_aa_seq(fields[8], qAaSeqLen);
        std::string t_aa_seq(fields[9], tAaSeqLen);
        std::string q_12st_seq(fields[10], q12stSeqLen);
        std::string t_12st_seq(fields[11], t12stSeqLen);

        // Parse score
        float score = atof(scoreStr.c_str());

        // Encode sequences from strings to numeric arrays
        unsigned char* q_aa_encoded = encodeSequence(q_aa_seq, ALPHABET_AA);
        unsigned char* t_aa_encoded = encodeSequence(t_aa_seq, ALPHABET_AA);
        unsigned char* q_3di_encoded = encodeSequence(q_3di_seq, ALPHABET_AA);  // 3Di uses same alphabet as AA
        unsigned char* t_3di_encoded = encodeSequence(t_3di_seq, ALPHABET_AA);
        unsigned char* q_12st_encoded = encodeSequence(q_12st_seq, ALPHABET_12ST);
        unsigned char* t_12st_encoded = encodeSequence(t_12st_seq, ALPHABET_12ST);

        // Initialize aligner with query sequences
        aligner.init_ssw(
            q_aa_encoded,
            q_3di_encoded,
            q_12st_encoded,
            q_aa_seq.length(),
            submat_aa,
            submat_3di,
            submat_12st,
            10000
        );

        // Perform alignment (gap_open=10, gap_extend=1, with backtrace)
        int computed_score = aligner.align(
            t_aa_encoded,
            t_3di_encoded,
            t_12st_encoded,
            t_aa_seq.length(),
            10,   // gap_open
            1,    // gap_extend
            true  // compute_backtrace
        );

        // Get backtrace and convert to CIGAR
        std::string backtrace = aligner.getBacktrace();
        std::string computed_cigar = backtraceToCigar(backtrace);

        // Get alignment start positions (0-based from aligner)
        int computed_qstart = aligner.getAlignStartRow();  // Keep as 0-based
        int computed_tstart = aligner.getAlignStartCol();  // Keep as 0-based

        // Parse expected start positions (0-based in test data)
        int expected_qstart = atoi(qStartStr.c_str());
        int expected_tstart = atoi(tStartStr.c_str());

        // Compare score, CIGAR, and start positions with expected values
        int score_int = (int)round(score);
        int score_diff = abs(computed_score - score_int);
        bool score_match = (score_diff == 0);
        bool cigar_match = (computed_cigar == expectedCigar);
        bool qstart_match = (computed_qstart == expected_qstart);
        bool tstart_match = (computed_tstart == expected_tstart);

        // Print test result (only show failures)
        if (score_match && cigar_match && qstart_match && tstart_match) {
            passedTests++;
        } else {
            std::cout << "Test " << totalTests << ": " << qname << " vs " << tname << " - FAIL" << std::endl;
            if (!score_match) {
                std::cout << "  Score: " << computed_score << " (expected: " << score_int << ", diff: " << score_diff << ")" << std::endl;
            } else {
                std::cout << "  Score: " << computed_score << " ✓" << std::endl;
            }
            if (!cigar_match) {
                std::cout << "  CIGAR expected: \"" << expectedCigar << "\"" << std::endl;
                std::cout << "  CIGAR computed: \"" << computed_cigar << "\"" << std::endl;
            } else {
                std::cout << "  CIGAR: " << computed_cigar << " ✓" << std::endl;
            }
            if (!qstart_match || !tstart_match) {
                std::cout << "  Start expected: q=" << expected_qstart << ", t=" << expected_tstart << std::endl;
                std::cout << "  Start computed: q=" << computed_qstart << ", t=" << computed_tstart << std::endl;
            } else {
                std::cout << "  Start: q=" << computed_qstart << ", t=" << computed_tstart << " ✓" << std::endl;
            }

            // Print aligned 3Di sequences for debugging (full length, no truncation)
            std::string q_aligned, t_aligned;
            getAlignedSequences(backtrace, q_3di_seq, t_3di_seq, computed_qstart, computed_tstart, q_aligned, t_aligned);
            std::cout << "  3Di alignment (length=" << q_aligned.length() << "):" << std::endl;
            std::cout << "    Q: " << q_aligned << std::endl;
            std::cout << "    T: " << t_aligned << std::endl;

            failedTests++;
        }

        // Cleanup sequences
        delete[] q_aa_encoded;
        delete[] t_aa_encoded;
        delete[] q_3di_encoded;
        delete[] t_3di_encoded;
        delete[] q_12st_encoded;
        delete[] t_12st_encoded;
    }

    fclose(fp);

    // Cleanup substitution matrices
    freeSubmat<21>(submat_aa);
    freeSubmat<21>(submat_3di);
    freeSubmat<13>(submat_12st);

    // Print summary
    std::cout << std::endl << "=== Test Summary ===" << std::endl;
    std::cout << "Total tests:  " << totalTests << std::endl;
    std::cout << "Passed:       " << passedTests << std::endl;
    std::cout << "Failed:       " << failedTests << std::endl;
    std::cout << "Success rate: " << (totalTests > 0 ? (100.0 * passedTests / totalTests) : 0.0) << "%" << std::endl;

    return (failedTests == 0) ? 0 : 1;
}
