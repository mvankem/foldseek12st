#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "StructureUtil.h"
#include "DistanceCalculator.h"
#include "QueryMatcher.h"
#include "Sequence.h"
#include "FastSort.h"

#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

// 12st ugaln submat (permutated)
static constexpr short PREFILTER_MATRIX[13][13] = {
    {  8,  3, -2,  2, -2, -4, -6, -8, -2, -5, -7, -5, -1 }, // A
    {  3,  6,  0, -2, -4, -3, -3, -5, -5, -3, -5, -7, -1 }, // C
    { -2,  0,  7,  4, -1,  1,  0, -3, -4, -2, -3, -5,  0 }, // D
    {  2, -2,  4,  8,  3, -2, -4, -7, -1, -3, -5, -2, -1 }, // E
    { -2, -4, -1,  3,  6,  2, -4, -6,  0, -2, -3, -2, -3 }, // F
    { -4, -3,  1, -2,  2,  6,  0, -3, -3, -1, -2, -4, -3 }, // G
    { -6, -3,  0, -4, -4,  0,  3, -2, -4,  0, -1, -5, -2 }, // H
    { -8, -5, -3, -7, -6, -3, -2,  0, -6, -2, -3, -8, -3 }, // I
    { -2, -5, -4, -1,  0, -3, -4, -6,  6,  2, -2,  1, -3 }, // K
    { -5, -3, -2, -3, -2, -1,  0, -2,  2,  5,  0, -2, -2 }, // L
    { -7, -5, -3, -5, -3, -2, -1, -3, -2,  0,  6,  3,  0 }, // M
    { -5, -7, -5, -2, -2, -4, -5, -8,  1, -2,  3,  7, -2 }, // N
    { -1, -1,  0, -1, -3, -3, -2, -3, -3, -2,  0, -2,  9 }  // X
};

// Compute ungapped alignment score for 12st alphabet only
template<typename T>
static int ungappedAlignment12st(const T *seq12st1,
                                  const T *seq12st2,
                                  const unsigned int length,
                                  short **sub12stMat) {
    int maxScore = 0;
    int score = 0;
    for(unsigned int pos = 0; pos < length; pos++){
        int curr12st = sub12stMat[static_cast<int>(seq12st1[pos])][static_cast<int>(seq12st2[pos])];
        score = curr12st + score;
        score = (score <= 0) ? 0 : score;
        maxScore = (score > maxScore) ? score : maxScore;
    }
    return maxScore;
}

int structurerescorediagonal12st(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    // Load 3Di alphabet for query and target
    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *t12stDbr = NULL;
    IndexReader *q12stDbr = NULL;

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        q12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        t12stDbr = q12stDbr;
    } else {
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        t12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        q12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    }

    if (q12stDbr == NULL || q12stDbr->sequenceReader == NULL) {
        Debug(Debug::ERROR) << "Failed to load query 12st database\n";
        EXIT(EXIT_FAILURE);
    }
    if (t12stDbr == NULL || t12stDbr->sequenceReader == NULL) {
        Debug(Debug::ERROR) << "Failed to load target 12st database\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    int dbtype = resultReader.getDbtype(); // Preserve input dbtype (DBTYPE_PREFILTER_RES)

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbtype);
    dbw.open();

    // TODO: The 12st substitution matrix is currently hardcoded (PREFILTER_MATRIX).
    // Consider loading from mat12st.out file for flexibility, similar to how 3Di matrix is loaded.

    // Build a short** view into the compile-time PREFILTER_MATRIX (no heap alloc)
    // PREFILTER_MATRIX is constexpr (const), so cast away const since we only read from it.
    static short *subMatrix12st_arr[13];
    short **subMatrix12st = subMatrix12st_arr;
    for (int i = 0; i < 13; ++i) {
        subMatrix12st[i] = const_cast<short*>(PREFILTER_MATRIX[i]);
    }

    Debug(Debug::INFO) << "Starting rescorediagonal12st for " << resultReader.getSize() << " queries\n";
    Debug(Debug::INFO) << "Using 12st substitution matrix (13x13: A-M + X)\n";

    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        // NOTE: We use 3Di substitution matrix for sequence mapping since the Sequence class
        // requires a BaseMatrix for character-to-index mapping. The actual scoring uses PREFILTER_MATRIX.
        SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);

        // Create Sequence objects for 12st alphabet mapping (similar to template)
        Sequence qSeq12st(par.maxSeqLen, q12stDbr->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, false);
        Sequence tSeq12st(par.maxSeqLen, t12stDbr->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, false);

        char buffer[1024+32768];
        std::string resultBuffer;
        std::vector<hit_t> hits;
        hits.reserve(300);

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);
            progress.updateProgress();

            if(*data != '\0') {
                unsigned int queryId = q12stDbr->sequenceReader->getId(queryKey);

                // Get 12st query sequence
                char *querySeq12st_raw = q12stDbr->sequenceReader->getData(queryId, thread_idx);
                unsigned int querySeqLen = q12stDbr->sequenceReader->getSeqLen(queryId);                

                // Map raw sequence to numeric representation (A,C,D,E -> 0,1,2,3)
                qSeq12st.mapSequence(queryId, queryKey, querySeq12st_raw, querySeqLen);
                while (*data != '\0') {
                    // Parse prefilter hit
                    hit_t prefHit = QueryMatcher::parsePrefilterHit(data);
                    data = Util::skipLine(data);

                    const unsigned int dbKey = prefHit.seqId;
                    unsigned int targetId = t12stDbr->sequenceReader->getId(dbKey);

                    // Get 12st target sequence
                    char *targetSeq12st_raw = t12stDbr->sequenceReader->getData(targetId, thread_idx);
                    const unsigned int targetSeqLen = t12stDbr->sequenceReader->getSeqLen(targetId);

                    // Map raw target sequence to numeric representation
                    tSeq12st.mapSequence(targetId, dbKey, targetSeq12st_raw, targetSeqLen);

                    // Compute diagonal alignment for 12st
                    // Note: prefHit.diagonal is unsigned short, but represents signed values via two's complement
                    int diagonal = static_cast<short>(prefHit.diagonal);
                    int score12st = 0;

                    unsigned int minDistToDiagonal = abs(diagonal);
                    if (diagonal >= 0 && minDistToDiagonal < querySeqLen) {
                        unsigned int minSeqLen = std::min(targetSeqLen, querySeqLen - minDistToDiagonal);
                        // Use numeric sequence (0,1,2,3...) and custom-loaded matrix
                        score12st = ungappedAlignment12st(qSeq12st.numSequence + minDistToDiagonal,
                                                          tSeq12st.numSequence,
                                                          minSeqLen, subMatrix12st);
                    } else if (diagonal < 0 && minDistToDiagonal < targetSeqLen) {
                        unsigned int minSeqLen = std::min(targetSeqLen - minDistToDiagonal, querySeqLen);
                        // Use numeric sequence (0,1,2,3...) and custom-loaded matrix
                        score12st = ungappedAlignment12st(qSeq12st.numSequence,
                                                          tSeq12st.numSequence + minDistToDiagonal,
                                                          minSeqLen, subMatrix12st);
                    }

                    // Add 12st score to existing 3Di score
                    int combinedScore = prefHit.prefScore + score12st;

                    // Collect hit for sorting
                    hit_t hit;
                    hit.seqId = dbKey;
                    hit.prefScore = combinedScore;
                    hit.diagonal = diagonal;
                    hits.emplace_back(hit);
                }
            }

            // Sort hits by score (descending) and seqId (ascending)
            if (hits.size() > 1) {
                SORT_SERIAL(hits.begin(), hits.end(), hit_t::compareHitsByScoreAndId);
            }

            // Write sorted hits to buffer
            for (size_t i = 0; i < hits.size(); ++i) {
                size_t len = QueryMatcher::prefilterHitToBuffer(buffer, hits[i]);
                resultBuffer.append(buffer, len);
            }

            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            hits.clear();
        }
    }

    dbw.close();
    resultReader.close();

    // No cleanup needed for static matrix

    if (sameDB == false) {
        delete t3DiDbr;
        delete t12stDbr;
    }
    delete q12stDbr;

    return EXIT_SUCCESS;
}
