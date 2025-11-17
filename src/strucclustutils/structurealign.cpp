#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "LDDT.h"
#include "Simple3AlphSW.h"

#ifdef OPENMP
#include <omp.h>
#endif

// need for sorting the results
static bool compareHitsByStructureBits(const Matcher::result_t &first, const Matcher::result_t &second) {
    if (first.score != second.score) {
        return first.score > second.score;
    }
    if (first.dbLen != second.dbLen) {
        return first.dbLen < second.dbLen;
    }
    return first.dbKey < second.dbKey;
}


static void structureAlignDefault(LocalParameters & par) {
    par.compBiasCorrectionScale = 0.5;
    par.alignmentType = LocalParameters::ALIGNMENT_TYPE_3DI_AA;
}

int alignStructure(StructureSmithWaterman & structureSmithWaterman,
                   StructureSmithWaterman & reverseStructureSmithWaterman,
                   Sequence & tSeqAA, Sequence & tSeq3Di,
                   unsigned int querySeqLen, unsigned int targetSeqLen,
                   EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                   Matcher::result_t & res, std::string & backtrace,
                   Parameters & par) {

    float seqId = 0.0;
    backtrace.clear();
    // align only score and end pos
    StructureSmithWaterman::s_align align = structureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                                    par.gapExtend.values.aminoacid(), querySeqLen / 2);
    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, align.qCov, align.tCov));
    if(hasLowerCoverage){
        return -1;
    }
    // we can already stop if this e-value isn't good enough, it wont be any better in the next step
    align.evalue = evaluer.computeEvalueCorr(align.score1, muLambda.first, muLambda.second);
    bool hasLowerEvalue = align.evalue > par.evalThr;
    if(hasLowerEvalue){
        return -1;
    }

    StructureSmithWaterman::s_align revAlign;
    //if(structureSmithWaterman.isProfileSearch()){
    //    revAlign.score1 = 0;
    //} else {
    revAlign = reverseStructureSmithWaterman.alignScoreEndPos<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence, tSeq3Di.numSequence,
                                                                  targetSeqLen, par.gapOpen.values.aminoacid(),
                                                                  par.gapExtend.values.aminoacid(), querySeqLen / 2);
    //}
    int32_t score = static_cast<int32_t>(align.score1) - static_cast<int32_t>(revAlign.score1);
    align.evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);
    hasLowerEvalue = align.evalue > par.evalThr;
    if (hasLowerEvalue) {
        return -1;
    }

    bool blockAlignFailed = false;
    if (structureSmithWaterman.isProfileSearch() == false) {
        StructureSmithWaterman::s_align alignTmp = structureSmithWaterman.alignStartPosBacktraceBlock(
            tSeqAA.numSequence, tSeq3Di.numSequence, targetSeqLen, par.gapOpen.values.aminoacid(),
            par.gapExtend.values.aminoacid(), backtrace, align
        );

        if (align.score1 == UINT32_MAX) {
            Debug(Debug::WARNING) << "block-align failed, falling back to normal alignment\n";
            blockAlignFailed = true;
        } else {
            align = alignTmp;
        }
    }

    if (blockAlignFailed || structureSmithWaterman.isProfileSearch()) {
        align = structureSmithWaterman.alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(tSeqAA.numSequence,
                                                                                               tSeq3Di.numSequence,
                                                                                               targetSeqLen,
                                                                                               par.gapOpen.values.aminoacid(),
                                                                                               par.gapExtend.values.aminoacid(),
                                                                                               par.alignmentMode,
                                                                                               backtrace, align,
                                                                                               par.covMode, par.covThr,
                                                                                               querySeqLen / 2);
    }

    unsigned int alnLength = Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
    if(backtrace.size() > 0){
        alnLength = backtrace.size();
        seqId = Util::computeSeqId(par.seqIdMode, align.identicalAACnt, querySeqLen, targetSeqLen, alnLength);
    }
    align.score1 = score;
    res = Matcher::result_t(tSeqAA.getDbKey(), align.score1, align.qCov, align.tCov, seqId, align.evalue, alnLength,
                            align.qStartPos1, align.qEndPos1, querySeqLen, align.dbStartPos1, align.dbEndPos1, targetSeqLen, backtrace);
    return 0;
}


int computeAlternativeAlignment(StructureSmithWaterman & structureSmithWaterman,
                                StructureSmithWaterman & reverseStructureSmithWaterman,
                                Sequence & tSeqAA, Sequence & tSeq3Di,
                                unsigned int querySeqLen, unsigned int targetSeqLen,
                                EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                                Matcher::result_t & result, Matcher::result_t & altRes,
                                std::string & backtrace, Parameters & par) {
    const unsigned char xAAIndex = tSeqAA.subMat->aa2num[static_cast<int>('X')];
    const unsigned char x3DiIndex = tSeq3Di.subMat->aa2num[static_cast<int>('X')];
    for (int pos = result.dbStartPos; pos < result.dbEndPos; ++pos) {
        tSeqAA.numSequence[pos] = xAAIndex;
        tSeq3Di.numSequence[pos] = x3DiIndex;
    }
    if (alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                       tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                       evaluer, muLambda, altRes, backtrace, par) == -1) {
        return -1;
    }
    if (Alignment::checkCriteria(altRes, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
        return 0;
    } else {
        return -1;
    }
}

// Helper function to create 12st substitution matrix (13x13)
// Returns int8_t** with matrix values scaled by the provided scale factor
int8_t** create12stSubMat(float scale) {
    // 12st matrix: 13x13 (ABCDEFGHIJKLX) - from TestSimple3AlphSW.cpp
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

    // Allocate int8_t matrix
    int8_t** discretized = new int8_t*[13];
    for (int i = 0; i < 13; i++) {
        discretized[i] = new int8_t[13];
        for (int j = 0; j < 13; j++) {
            // Apply scaling and round to int8_t
            float scaled = SUBMAT_12ST_FLOAT[i][j] * scale;
            discretized[i][j] = (int8_t)round(scaled);  // FOR DEBUGGING PURPOSES!
        }
    }
    return discretized;
}

// Helper function to create 2D pointer array from flat 1D array
// This allows using flat arrays (tinySubMatAA, tinySubMat3Di) as 2D matrices
int8_t** create2DMatrixFromFlat(int8_t* flatMatrix, int size) {
    int8_t** matrix2D = new int8_t*[size];
    for (int i = 0; i < size; i++) {
        matrix2D[i] = &flatMatrix[i * size];
    }
    return matrix2D;
}

int alignStructure12st(Simple3AlphSW & simple3AlphSW,
                       Simple3AlphSW & reverseSimple3AlphSW,
                       Sequence & tSeqAA, Sequence & tSeq3Di, Sequence & tSeq12st,
                       unsigned int querySeqLen, unsigned int targetSeqLen,
                       EvalueNeuralNet & evaluer, std::pair<double, double> muLambda,
                       Matcher::result_t & res, std::string & backtrace,
                       Parameters & par) {
    float seqId = 0.0;
    backtrace.clear();

    // Perform alignment (backtrace is computed from DP matrices, no need for compute_backtrace flag)
    int score = simple3AlphSW.align(
        tSeqAA.numSequence,
        tSeq3Di.numSequence,
        tSeq12st.numSequence,
        targetSeqLen,
        par.gapOpen.values.aminoacid(),
        par.gapExtend.values.aminoacid(),
        false  // compute_backtrace not needed - backtrace() recomputes from DP matrices
    );

    // Get alignment end positions from max score position
    // maxScoreRow/Col are in DP matrix coordinates (1-indexed)
    // Convert to 0-based sequence coordinates
    unsigned int qEndPos = simple3AlphSW.getMaxScoreRow() - 1;
    unsigned int dbEndPos = simple3AlphSW.getMaxScoreCol() - 1;

    // Calculate coverage based on end positions
    // For initial check, assume alignment starts at position 0
    float qCov = static_cast<float>(qEndPos + 1) / static_cast<float>(querySeqLen);
    float tCov = static_cast<float>(dbEndPos + 1) / static_cast<float>(targetSeqLen);

    // Check coverage
    bool hasLowerCoverage = !(Util::hasCoverage(par.covThr, par.covMode, qCov, tCov));
    if(hasLowerCoverage){
        return -1;
    }

    // Compute E-value with initial score
    // TODO: need to rescor alignment only with 3Di + AA
    double evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);
    evalue = 0;  // change this
    bool hasLowerEvalue = evalue > par.evalThr;
    if(hasLowerEvalue){
        return -1;
    }

    // Reverse alignment for score adjustment (like alignStructure lines 62-75)
    int reverseScore = reverseSimple3AlphSW.align(
        tSeqAA.numSequence,
        tSeq3Di.numSequence,
        tSeq12st.numSequence,
        targetSeqLen,
        par.gapOpen.values.aminoacid(),
        par.gapExtend.values.aminoacid(),
        false  // No backtrace needed for reverse alignment
    );

    // Adjust score: forward - reverse (like alignStructure line 70)
    int adjustedScore = score - reverseScore;

    // Recompute E-value with adjusted score (like alignStructure line 71)
    evalue = evaluer.computeEvalueCorr(adjustedScore, muLambda.first, muLambda.second);
    // TODO: need to rescor alignment only with 3Di + AA
    evalue = 0;  // change this

    // Check E-value with adjusted score (like alignStructure lines 72-75)
    hasLowerEvalue = evalue > par.evalThr;
    if(hasLowerEvalue){
        return -1;
    }

    // Use adjusted score for final result (like alignStructure line 109)
    score = adjustedScore;

    // Get backtrace string (M/I/D operations) - computed from DP matrices
    backtrace = simple3AlphSW.getBacktrace();

    // Get alignment start positions (0-based) - computed during backtrace
    unsigned int qStartPos = simple3AlphSW.getAlignStartRow();
    unsigned int dbStartPos = simple3AlphSW.getAlignStartCol();

    // Recalculate coverage based on actual alignment span
    qCov = static_cast<float>(qEndPos - qStartPos + 1) / static_cast<float>(querySeqLen);
    tCov = static_cast<float>(dbEndPos - dbStartPos + 1) / static_cast<float>(targetSeqLen);

    // Calculate alignment length from backtrace
    unsigned int alnLength = backtrace.size();

    // Calculate sequence identity
    // TODO: Count identical AA positions in alignment
    unsigned int identicalAACnt = alnLength; // Placeholder, should be computed properly
    if(backtrace.size() > 0){
        seqId = Util::computeSeqId(par.seqIdMode, identicalAACnt, querySeqLen, targetSeqLen, alnLength);
    }

    // Construct Matcher::result_t with all alignment details
    res = Matcher::result_t(tSeqAA.getDbKey(), score, qCov, tCov, seqId, evalue, alnLength,
                            qStartPos, qEndPos, querySeqLen, dbStartPos, dbEndPos, targetSeqLen, backtrace);
    return 0;
}


int structurealign(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    structureAlignDefault(par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool sameDB = false;
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool alignmentIsExtended = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;
    IndexReader tAADbr(par.db2, par.threads,
                             alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                             (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()),
                                              Parameters::DBTYPE_INDEX_DB);

    IndexReader t3DiDbr(is3DiIdx ? t3DiDbrName : par.db2, par.threads,
                              alignmentIsExtended ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                              (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                              DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                              alignmentIsExtended ? "_seq_ss" : "_ss");

    IndexReader *q3DiDbr = NULL;
    IndexReader *qAADbr = NULL;

    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        q3DiDbr = &t3DiDbr;
        qAADbr = &tAADbr;
    } else {
        qAADbr = new IndexReader(par.db1, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
        q3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    // Load 12st sequences for alignment type 4 (3Di+AA+12st)
    IndexReader *q12stDbr = NULL;
    IndexReader *t12stDbr = NULL;
    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
        if (sameDB) {
            q12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
            t12stDbr = q12stDbr;
        } else {
            q12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db1, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
            t12stDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss2"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        }
    }

    bool db1CaExist = FileUtil::fileExists((par.db1 + "_ca.dbtype").c_str());
    bool db2CaExist = FileUtil::fileExists((par.db2 + "_ca.dbtype").c_str());
    if(Parameters::isEqualDbtype(tAADbr.getDbtype(), Parameters::DBTYPE_INDEX_DB)){
        db2CaExist = true;
    }
    if(par.sortByStructureBits) {
        bool disableStructureBits = false;
        if(db1CaExist == false || db2CaExist == false){
            Debug(Debug::WARNING) << "Cannot find " << FileUtil::baseName(par.db1) << " C-alpha or " << FileUtil::baseName(par.db2) << " C-alpha database\n";
            disableStructureBits = true;
        }
        if(par.alignmentMode == 1 || par.alignmentMode == 2){
            Debug(Debug::WARNING) << "Cannot use --sort-by-structure-bits 1 with --alignment-mode 1 or 2\n";
            disableStructureBits = true;
        }
        if(disableStructureBits){
            Debug(Debug::WARNING) << "Disabling --sort-by-structure-bits\n";
            Debug(Debug::WARNING) << "This impacts the final score and ranking of hits, but not E-values themselves. Ranking alterations primarily occur for E-values < 10^-1.\n";
            par.sortByStructureBits = false;
        }
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbtype =  Parameters::DBTYPE_ALIGNMENT_RES;
    if(alignmentIsExtended){
        dbtype = DBReader<unsigned int>::setExtendedDbtype(dbtype, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  dbtype);
    dbw.open();

    bool needTMaligner = (par.tmScoreThr > 0);
    bool needLDDT = (par.lddtThr > 0);
    if (par.sortByStructureBits) {
        needLDDT = true;
        needTMaligner = true;
    } else {
        if (needTMaligner && (db1CaExist == false || db2CaExist == false)) {
            Debug(Debug::WARNING) << "Cannot use --tmscore-threshold with --sort-by-structure-bits 0\n"
                                  << "Disabling --tmscore-threshold\n";
            needTMaligner = false;
        }
        if (needLDDT && (db1CaExist == false || db2CaExist == false)) {
            Debug(Debug::WARNING) << "Cannot use --lddt-threshold with --sort-by-structure-bits 0\n"
                                  << "Disabling --lddt-threshold\n";
            needLDDT = false;
        }
    }
    bool needCalpha = (needTMaligner || needLDDT);
    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;
    if(needCalpha){
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
             tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    alignmentIsExtended ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) :
                                           IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                    touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    alignmentIsExtended ? "_seq_ca" : "_ca"
            );
        }
    }

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    float aaFactor = (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_AA ||
                      par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) ? 1.4 : 0.0;
    SubstitutionMatrix subMatAA(blosum.c_str(), aaFactor, par.scoreBias);
    //temporary output file
    Debug::Progress progress(resultReader.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);

    for (int i = 0; i < subMat3Di.alphabetSize; i++) {
        for (int j = 0; j < subMat3Di.alphabetSize; j++) {
            tinySubMat3Di[i * subMat3Di.alphabetSize + j] = subMat3Di.subMatrix[i][j]; // for farrar profile
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }

    // Initialize 12st substitution matrix if alignment type is 4
    int8_t** submat_12st = NULL;
    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
        submat_12st = create12stSubMat(par.submat12stScale);
    }

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr.sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale, &subMatAA, &subMat3Di);
        Simple3AlphSW simple3AlphSW;  // For alignment type 4 (3Di+AA+12st)
        Simple3AlphSW reverseSimple3AlphSW;  // For reverse alignment scoring
        // Create 2D matrices for Simple3AlphSW from flat arrays
        int8_t** tinySubMatAA_2d = NULL;
        int8_t** tinySubMat3Di_2d = NULL;
        if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
            tinySubMatAA_2d = create2DMatrixFromFlat(tinySubMatAA, subMatAA.alphabetSize);
            tinySubMat3Di_2d = create2DMatrixFromFlat(tinySubMat3Di, subMat3Di.alphabetSize);
        }
        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(q3DiDbr->sequenceReader->getMaxSeqLen() + 1, t3DiDbr.sequenceReader->getMaxSeqLen() + 1), false, true, par.exactTMscore);
        }
        LDDTCalculator *lddtcalculator = NULL;
        if(needLDDT) {
            lddtcalculator = new LDDTCalculator(q3DiDbr->sequenceReader->getMaxSeqLen() + 1,  t3DiDbr.sequenceReader->getMaxSeqLen() + 1);
        }
        Sequence qSeqAA(par.maxSeqLen, qAADbr->getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, q3DiDbr->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        // 12st sequences for alignment type 4 (using 3Di submat temporarily)
        Sequence *qSeq12st = NULL;
        Sequence *tSeq12st = NULL;
        if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
            qSeq12st = new Sequence(par.maxSeqLen, q12stDbr->getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
            tSeq12st = new Sequence(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        }
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        TMaligner::TMscoreResult tmres;
        LDDTCalculator::LDDTScoreResult lddtres;
        // write output file

#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();
            char *data = resultReader.getData(id, thread_idx);
            size_t queryKey = resultReader.getDbKey(id);
            if(*data != '\0') {
                unsigned int queryId = q3DiDbr->sequenceReader->getId(queryKey);

                char *querySeqAA = qAADbr->sequenceReader->getData(queryId, thread_idx);
                char *querySeq3Di = q3DiDbr->sequenceReader->getData(queryId, thread_idx);
                unsigned int querySeqLen = q3DiDbr->sequenceReader->getSeqLen(queryId);
                qSeq3Di.mapSequence(id, queryKey, querySeq3Di, querySeqLen);
                qSeqAA.mapSequence(id, queryKey, querySeqAA, querySeqLen);
                // Map 12st sequences for alignment type 4
                if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
                    unsigned int q12stId = q12stDbr->sequenceReader->getId(queryKey);
                    char *querySeq12st = q12stDbr->sequenceReader->getData(q12stId, thread_idx);
                    unsigned int querySeqLen12st = q12stDbr->sequenceReader->getSeqLen(q12stId);
                    qSeq12st->mapSequence(id, queryKey, querySeq12st, querySeqLen12st);
                }
                if(needCalpha){
                    size_t qId = qcadbr->sequenceReader->getId(queryKey);
                    char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
                    size_t qCaLength = qcadbr->sequenceReader->getEntryLen(qId);
                    float* queryCaData = qcoords.read(qcadata, qSeq3Di.L, qCaLength);
                    if(needTMaligner){
                        tmaligner->initQuery(queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L], NULL, qSeq3Di.L);
                    }
                    if(needLDDT){
                        lddtcalculator->initQuery(qSeq3Di.L, queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L+qSeq3Di.L]);
                    }
                }
                std::pair<double, double> muLambda = evaluer.predictMuLambda(qSeq3Di.numSequence, qSeq3Di.L);
                // Skip StructureSmithWaterman initialization for alignment type 4 (3Di+AA+12st)
                if (par.alignmentType != LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
                    structureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                    qSeq3Di.reverse();
                    qSeqAA.reverse();
                    reverseStructureSmithWaterman.ssw_init(&qSeqAA, &qSeq3Di, tinySubMatAA, tinySubMat3Di, &subMatAA);
                } else {
                    // Initialize Simple3AlphSW for alignment type 4 (3Di+AA+12st)
                    // First initialize forward aligner
                    simple3AlphSW.init_ssw(
                        qSeqAA.numSequence,
                        qSeq3Di.numSequence,
                        qSeq12st->numSequence,
                        querySeqLen,
                        tinySubMatAA_2d,
                        tinySubMat3Di_2d,
                        submat_12st,
                        t3DiDbr.sequenceReader->getMaxSeqLen()
                    );
                    // Reverse query sequences for reverse aligner
                    qSeq3Di.reverse();
                    qSeqAA.reverse();
                    qSeq12st->reverse();
                    // Initialize reverse aligner with reversed sequences
                    reverseSimple3AlphSW.init_ssw(
                        qSeqAA.numSequence,
                        qSeq3Di.numSequence,
                        qSeq12st->numSequence,
                        querySeqLen,
                        tinySubMatAA_2d,
                        tinySubMat3Di_2d,
                        submat_12st,
                        t3DiDbr.sequenceReader->getMaxSeqLen()
                    );
                }
                int passedNum = 0;
                int rejected = 0;
                while (*data != '\0' && passedNum < par.maxAccept && rejected < par.maxRejected) {
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    data = Util::skipLine(data);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = t3DiDbr.sequenceReader->getId(dbKey);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    char * targetSeq3Di = t3DiDbr.sequenceReader->getData(targetId, thread_idx);
                    char * targetSeqAA = tAADbr.sequenceReader->getData(targetId, thread_idx);
                    const int targetSeqLen = static_cast<int>(t3DiDbr.sequenceReader->getSeqLen(targetId));

                    tSeq3Di.mapSequence(targetId, dbKey, targetSeq3Di, targetSeqLen);
                    tSeqAA.mapSequence(targetId, dbKey, targetSeqAA, targetSeqLen);
                    // Map 12st target sequences for alignment type 4
                    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
                        unsigned int t12stId = t12stDbr->sequenceReader->getId(dbKey);
                        char *targetSeq12st = t12stDbr->sequenceReader->getData(t12stId, thread_idx);
                        unsigned int targetSeqLen12st = t12stDbr->sequenceReader->getSeqLen(t12stId);
                        tSeq12st->mapSequence(targetId, dbKey, targetSeq12st, targetSeqLen12st);
                    }
                    if(Util::canBeCovered(par.covThr, par.covMode, qSeq3Di.L, targetSeqLen) == false){
                        rejected++;
                        continue;
                    }
                    Matcher::result_t res;
                    // Alignment for type 4 (3Di+AA+12st) uses Simple3AlphSW
                    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
                        if(alignStructure12st(simple3AlphSW, reverseSimple3AlphSW,
                                            tSeqAA, tSeq3Di, *tSeq12st,
                                            querySeqLen, targetSeqLen,
                                            evaluer, muLambda, res, backtrace, par) == -1){
                            rejected++;
                            continue;
                        }
                    } else {
                        if(alignStructure(structureSmithWaterman, reverseStructureSmithWaterman,
                                        tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                        evaluer, muLambda, res, backtrace, par) == -1){
                            rejected++;
                            continue;
                        }
                    }

                    if (Alignment::checkCriteria(res, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        if(needCalpha) {
                            size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                            char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                            size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                            float* targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                            if(needTMaligner) {
                                tmres = tmaligner->computeTMscore(targetCaData,
                                                                  &targetCaData[res.dbLen],
                                                                  &targetCaData[res.dbLen +
                                                                                res.dbLen],
                                                                  res.dbLen,
                                                                  res.qStartPos,
                                                                  res.dbStartPos,
                                                                  res.backtrace,
                                                                  TMaligner::normalization(par.tmScoreThrMode, std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos ), res.qLen, res.dbLen));

                                if (tmres.tmscore < par.tmScoreThr) {
                                    continue;
                                }
                            }
                            if(needLDDT){
                                lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos,
                                                                           res.backtrace,
                                                                           targetCaData, &targetCaData[res.dbLen],
                                                                           &targetCaData[res.dbLen+res.dbLen]);

                                if(lddtres.avgLddtScore < par.lddtThr){
                                    continue;
                                }
                                res.dbcov = lddtres.avgLddtScore;
                            }
                            if(par.sortByStructureBits && needTMaligner && needLDDT){
                                res.score = res.score * sqrt(lddtres.avgLddtScore * tmres.tmscore);
                            }
                        }


                        alignmentResult.emplace_back(res);
                        int altAli = par.altAlignment;
                        bool moreAltAli = true;
                        while(altAli && moreAltAli){
                            Matcher::result_t altRes;
                            if(computeAlternativeAlignment(structureSmithWaterman, reverseStructureSmithWaterman,
                                                           tSeqAA, tSeq3Di, querySeqLen, targetSeqLen,
                                                           evaluer, muLambda, res, altRes,
                                                           backtrace, par) == -1) {
                                moreAltAli = false;
                                continue;
                            }
                            alignmentResult.push_back(altRes);
                            res = altRes;
                            altAli--;
                        }
                        passedNum++;
                        rejected = 0;
                    } else {
                        rejected++;
                    }
                }
            }


            if (alignmentResult.size() > 1) {
                if(par.sortByStructureBits) {
                    SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), compareHitsByStructureBits);
                } else {
                    SORT_SERIAL(alignmentResult.begin(), alignmentResult.end(), Matcher::compareHits);
                }
            }
            for (size_t result = 0; result < alignmentResult.size(); result++) {
                size_t len = Matcher::resultToBuffer(buffer, alignmentResult[result], par.addBacktrace);
                resultBuffer.append(buffer, len);
            }
            dbw.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
            resultBuffer.clear();
            alignmentResult.clear();
        }
        if(needTMaligner){
            delete tmaligner;
        }
        if(needLDDT){
            delete lddtcalculator;
        }
        if(par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA){
            delete qSeq12st;
            delete tSeq12st;
        }
        // Cleanup 2D matrices (only pointer arrays, not the underlying data)
        if (tinySubMatAA_2d != NULL) {
            delete[] tinySubMatAA_2d;
        }
        if (tinySubMat3Di_2d != NULL) {
            delete[] tinySubMat3Di_2d;
        }
    }

    free(tinySubMatAA);
    free(tinySubMat3Di);

    // Cleanup 12st substitution matrix if it was allocated
    if (submat_12st != NULL) {
        for (int i = 0; i < 13; i++) {
            delete[] submat_12st[i];
        }
        delete[] submat_12st;
    }

    dbw.close();
    resultReader.close();

    if(needCalpha){
        if (sameDB == false) {
            delete tcadbr;
        }
        delete qcadbr;
    }

    if (sameDB == false) {
        delete q3DiDbr;
        delete qAADbr;
    }

    // Cleanup 12st IndexReaders for alignment type 4
    if (par.alignmentType == LocalParameters::ALIGNMENT_TYPE_3DI_12ST_AA) {
        if (sameDB == false) {
            delete t12stDbr;
        }
        delete q12stDbr;
    }

    return EXIT_SUCCESS;
}
