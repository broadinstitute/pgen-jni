#define BOOST_TEST_MODULE pgen_write
#include <sys/stat.h>
#include <stdio.h>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/array.hpp>
#include "pgenException.h"
#include "pgenMissingVariantsException.h"
#include "pgenEmptyPgenException.h"
#include "pgenContext.h"
#include "pgenIO.h"
#include "pgenUtils.h"

using namespace boost::unit_test;
using namespace pgenlib;

// Unit level tests for the PGEN writer. Very little validation of the resulting pgen files is done here; the
// Java tests in the enclosing Java project do the actual validation and round trip concordance verification.

//******************* Forward Declarations/Constants *******************
constexpr int TMP_FILENAME_SIZE = 4096;
template<size_t N> void CreateTempFile(const char* const nameTemplate, char (&outputFileName)[N]);
void GenerateAlleleCodeDistribution(int32_t* const allele_codes, const int n_samples, const int n_alleles);
long WriteTestPgen(
        const int32_t* const allele_codes,
        const unsigned char* const phase_bytes,
        const int32_t allele_ct,
        const uint32_t pgen_file_mode,
        const uint32_t write_flags,
        const long n_variants,
        const int n_samples,
        long &writtenVariantCount);
// integer constants to parallel PgenFileMode, for use when calling jni callable functions, which can't
// use the PgenFileMode enum provided by plink2
constexpr uint32_t PGEN_FILE_MODE_BACKWARD_SEEK = static_cast<int>(plink2::PgenWriteMode::kPgenWriteBackwardSeek);
constexpr uint32_t PGEN_FILE_MODE_WRITE_SEPARATE_INDEX = static_cast<int>(plink2::PgenWriteMode::kPgenWriteSeparateIndex);
constexpr uint32_t PGEN_FILE_MODE_WRITE_AND_COPY = static_cast<int>(plink2::PgenWriteMode::kPgenWriteAndCopy);

//******************* Tests *******************

// simple test to exercise throwing/catching of PgenException
BOOST_AUTO_TEST_CASE(TestExceptionPropagation) {
    const char* const expectedPropagationMessage = "Fake pgen exception";
    BOOST_REQUIRE_EXCEPTION(
            throw PgenException(expectedPropagationMessage),
            PgenException,
            [expectedPropagationMessage](PgenException ex) -> bool {
                return !strcmp(ex.what(), expectedPropagationMessage);
            }
    );
}

// verify for at least one case that we can sccessfully convert a plink2::pglerr enum into a string
BOOST_AUTO_TEST_CASE(TestPglerrStringConversion) {
    const char* const expectedMessage = "kPglRetNotYetSupported";
    BOOST_REQUIRE_EXCEPTION(
            throwOnPglErr(plink2::PglErr::ec::kPglRetNotYetSupported, "Testing PglErr conversion"),
            PgenException,
            [expectedMessage](PgenException ex) -> bool {
                return strstr(ex.what(), expectedMessage);
            }
    );
}

// write a small, bi-allelic pgen file, once with each possible file write mode
constexpr boost::array<uint32_t, 3> s_pgenFileMode {
    PGEN_FILE_MODE_BACKWARD_SEEK,
    PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
    PGEN_FILE_MODE_WRITE_AND_COPY
};
BOOST_DATA_TEST_CASE(TestUnphasedBiallelicSmall, s_pgenFileMode) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples
    constexpr int32_t allele_codes[] {0, 0, 0, 1, 1, 1 };
    long variantCount = 0L;
    // don't be fooled by the reference to the variable "sample" below; its a variable name introduced by
    // the boost macro to refer to the parameter for the test case, which in this test is the pgen file mode...
    const long file_size = WriteTestPgen(allele_codes, nullptr, 2, sample, 0, n_variants, n_samples, variantCount);
    // only check for non-zero file size, since the file size varies with the file mode
    BOOST_REQUIRE_NE(file_size, 0);
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// write a larger, bi-allelic pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(TestUnphasedBiallelicLarge) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 2;
    // one variants's worth of allele codes drawn from 2 alleles
    int32_t *allele_codes = new int32_t[n_samples * 2];
    GenerateAlleleCodeDistribution(allele_codes, n_samples, n_alleles);
    long variantCount = 0L;
    const long file_size = WriteTestPgen(allele_codes, nullptr,n_alleles,PGEN_FILE_MODE_WRITE_AND_COPY, 0, n_variants, n_samples, variantCount);
    delete[] allele_codes;

    BOOST_REQUIRE_EQUAL(file_size, 125450028); // cause thats what it is
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// write a larger, multi-allelic, unphased pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(TestUnphasedMultiAllelicLarge) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 7;
    // synthesize one variants's worth of allele codes, with genotypes randomly drawn from 7 allele codes
    int32_t *allele_codes = new int32_t[n_samples * 2];
    GenerateAlleleCodeDistribution(allele_codes, n_samples, n_alleles);
    long variantCount = 0L;
    const long file_size = WriteTestPgen(allele_codes, nullptr, n_alleles,PGEN_FILE_MODE_WRITE_AND_COPY, 0, n_variants, n_samples, variantCount);
    delete[] allele_codes;

    BOOST_REQUIRE_EQUAL(file_size, 911252530); // cause thats what it is
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// write a larger, multi-allelic, phased pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(TestPhasedMultiAllelicLarge) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 7;
    // synthesize one variants's worth of allele codes, with genotypes randomly drawn from 7 allele codes
    int32_t *allele_codes = new int32_t[n_samples * 2];
    GenerateAlleleCodeDistribution(allele_codes, n_samples, n_alleles);
    unsigned char *phase_bytes = new unsigned char[n_samples];
    for (int i = 0; i < n_samples; i++) {
        phase_bytes[i] = (unsigned char) 0x1;
    }
    long variantCount = 0L;
    const long file_size = WriteTestPgen(
            allele_codes,
            phase_bytes,
            n_alleles,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            pgenlib::kWriteFlagMultiAllelic | pgenlib::kWriteFlagPreservePhasing,
            n_variants,
            n_samples,
            variantCount);
    delete[] allele_codes;
    delete[] phase_bytes;

    BOOST_REQUIRE_EQUAL(file_size, 1036402530); // cause thats what it is
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// verify that the issue described here is fixed: https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ for bi-allelics
BOOST_AUTO_TEST_CASE(TestBiallelicOneAlleleNotObserved) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    constexpr int n_alleles = 2;
    // one variants's worth of allele codes - drawn from 2 alleles over 3 samples, but with only one
    // allele actually observed
    constexpr int32_t allele_codes[] {0, 0, 0, 0, 0, 0 };
    long variantCount = 0L;
    WriteTestPgen(allele_codes, nullptr, n_alleles, PGEN_FILE_MODE_WRITE_AND_COPY, 0, n_variants, n_samples, variantCount);
}

// verify issue described here: https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ for multi-allelics
BOOST_AUTO_TEST_CASE(TestMultiallelicSomeAllelesNotObserved) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    constexpr int n_alleles = 7;
    // one variants's worth of allele codes - drawn from 7 alleles over 3 samples, but with only 3
    // alleles actually observed
    constexpr int32_t allele_codes[] {0, 1, 0, 5, 0, 4 };
    long variantCount = 0L;
    WriteTestPgen(allele_codes, nullptr, n_alleles, PGEN_FILE_MODE_WRITE_AND_COPY, 0, n_variants, n_samples, variantCount);
}

// claim that we're going to write 10 variants, but don't write any)
BOOST_AUTO_TEST_CASE(TestCloseNoWriteKnownVariantCount) {
    const char* const expectedNoWriteMessage = "closePgen called with number of variants written";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::OpenPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            false,
            10L,
            3,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            ClosePgen(pgenContext, 0),
            PgenMissingVariantsException,
            [expectedNoWriteMessage](PgenMissingVariantsException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
}

// unknown variant count, don't write any)
BOOST_AUTO_TEST_CASE(TestCloseNoWritesUnknownVariantCount) {
    const char* const expectedNoWriteMessage = "An empty PGEN is not valid";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::OpenPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            0,
            static_cast<long>(plink2::kPglMaxVariantCt),
            3,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            ClosePgen(pgenContext, 0),
            PgenEmptyPgenException,
            [expectedNoWriteMessage](PgenEmptyPgenException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
}

// claim that we're going to write 10 variants, but only write 1
BOOST_AUTO_TEST_CASE(TestCloseTooFewWritesKnownVariantCount) {
    constexpr int n_samples = 3;
    const char* const expectedNoWriteMessage = "closePgen called with number of variants written";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::OpenPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            0,
            10L,
            3,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    int32_t *allele_codes = new int32_t[n_samples * 2]{0};
    pgenlib::AppendAlleles(pgenContext, allele_codes, nullptr, 2);
    BOOST_REQUIRE_EXCEPTION(
            ClosePgen(pgenContext, 0),
            PgenMissingVariantsException,
            [expectedNoWriteMessage](PgenMissingVariantsException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
    delete[] allele_codes;
}

// claim that we're going to write 10 variants, but only write 1
BOOST_AUTO_TEST_CASE(TestCloseTooFewWritesUnknownVariantCount) {
    constexpr int n_samples = 3;
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::OpenPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            0,
            static_cast<long>(plink2::kPglMaxVariantCt),
            n_samples,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    int32_t *allele_codes = new int32_t[n_samples * 2]{0};
    pgenlib::AppendAlleles(pgenContext, allele_codes, nullptr, 2);
    ClosePgen(pgenContext, 0);
    delete[] allele_codes;
}

BOOST_AUTO_TEST_CASE(TestRejectInvalidAlleleCode) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples, with on bad allele code
    constexpr int32_t allele_codes[] {0, 0, 0, -17, 0, 0};
    const char* const expectedMessageFragment = "Attempt to append invalid allele code";
    long variantCount;
    BOOST_REQUIRE_EXCEPTION(
            WriteTestPgen(allele_codes, nullptr, 2, PGEN_FILE_MODE_WRITE_AND_COPY, 0, n_variants, n_samples, variantCount),
            PgenException,
            [expectedMessageFragment](PgenException ex) -> bool {
                return strstr(ex.what(), expectedMessageFragment);
            }
    );
}

// use a bogus pgen write mode
BOOST_AUTO_TEST_CASE(TestRejectInvalidWriteMode) {
    const char* const expectedInvalidModeMessage = "Invalid pgenWriteMode";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const int invalidWriteMode = -2; // must be one of 0, 1, 2
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::OpenPgen(tmpFileName, invalidWriteMode, 0, 10L, 3, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(TestRejectInvalidSampleCount) {
    const char* const expectedInvalidSampleCountMessage = "Invalid sample count";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const int invalidSampleCount = 0; // must be >= 1
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::OpenPgen(tmpFileName, PGEN_FILE_MODE_WRITE_AND_COPY, 0, 10L, invalidSampleCount, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidSampleCountMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidSampleCountMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(TestRejectInvalidVariantCount) {
    const char* const expectedInvalidVariantCountMessage = "Invalid variant count";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const long invalidVariantCount = 0; // must be >= 1
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::OpenPgen(tmpFileName, PGEN_FILE_MODE_WRITE_AND_COPY, false, invalidVariantCount, 3, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidVariantCountMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidVariantCountMessage);
            }
    );
}

// use pgen write mode 0 (seek) with an unknown variant count
BOOST_AUTO_TEST_CASE(TestRejectSeekWriteModeWithUnknownVariantCount) {
    const char* const expectedInvalidModeMessage = "requires a known variant count";
    char tmpFileName[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::OpenPgen(tmpFileName,
                              PGEN_FILE_MODE_BACKWARD_SEEK,
                              0,
                              static_cast<long>(plink2::kPglMaxVariantCt),
                              3,
                              plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

//******************* Local Test Utilities *******************

// generate allele codes for n_samples, distributed across values from n_alleles
void GenerateAlleleCodeDistribution(int32_t* const allele_codes, const int n_samples, const int n_alleles) {
    for (int i = 0; i < (n_samples * 2); i+=2) {
        allele_codes[i] = i % n_alleles;
        allele_codes[i+1] = (i + 1) % n_alleles;
    }
}

// the caller should call unlink() on the resulting file to cause it to be deleted
template<size_t N>
void CreateTempFile(const char* const nameTemplate, char (&outputFileName)[N]) {
    //this is deprecated (and maybe a little sketchy), but works nicely to obtain a tmp dir location
    std::string tmpPath = std::tmpnam(nullptr);
    snprintf(outputFileName, N, "%s_pgenBoostXXXXXX%s", tmpPath.c_str(), nameTemplate);

    // we don't actually need to create the file here, just reserve it
    int fDesc = mkstemps(outputFileName, strlen(nameTemplate));
    if (fDesc < 1) {
        char errMessage[pgenlib::kReservedMessageBufSize];
        snprintf(errMessage,
                 pgenlib::kReservedMessageBufSize,
                 "Temp file creation failed for (%s) with error(%s)",
                 outputFileName,
                 strerror(errno));
        throw PgenException(errMessage);
    }
    // the file is open, so close it, but leave the call to unlink so the caller can control when it is deleted
    close(fDesc);
}

// write a PGEN file given allele codes (the same allele code vector is used for each variant) and phase_bytes
// (may be null), # of variants, # of samples, and write mode
long WriteTestPgen(
        const int32_t *allele_codes,
        const unsigned char *phase_bytes,
        const int32_t allele_ct,
        const uint32_t pgen_file_mode,
        const uint32_t writeFlags,
        const long n_variants,
        const int n_samples,
        long &variantCount) {
    char tmp_file_name[TMP_FILENAME_SIZE];
    CreateTempFile("test_write.pgen", tmp_file_name);

    const pgenlib::PgenContext *const pgen_context = pgenlib::OpenPgen(
            tmp_file_name,
            pgen_file_mode,
            writeFlags,
            n_variants,
            n_samples,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgen_context, nullptr);

    for (int i = 0; i < n_variants; i++) {
        pgenlib::AppendAlleles(pgen_context, allele_codes, phase_bytes, allele_ct);
    }
    variantCount = GetNumberOfVariantsWritten(pgen_context);
    ClosePgen(pgen_context, 0);

    // for now, just validate that the file has SOME contents; the enclosing pgen project has test code
    // that verifies the contents using plink2
    struct stat st;
    stat(tmp_file_name, &st);
    const long file_size = st.st_size;
    BOOST_REQUIRE_NE(file_size, 0);

    unlink(tmp_file_name);
    return file_size;
}
