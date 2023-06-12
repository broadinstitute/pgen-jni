#define BOOST_TEST_MODULE pgen_write
#include <sys/stat.h>
#include <stdio.h>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/array.hpp>
#include "pgenException.h"
#include "pgenContext.h"
#include "pgenIO.h"
#include "pgenUtils.h"

using namespace boost::unit_test;
using namespace pgenlib;

// Unit level tests for the PGEN writer. Very little validation of the resulting pgen files is done here; the
// Java tests in the enclosing Java project do the actual validation and round trip concordance verification.

//******************* Forward Declarations/Constants *******************
constexpr int TMP_FILENAME_SIZE = 4096;
template<size_t N> void createTempFile(const char* const nameTemplate, char (&outputFileName)[N]);
void generate_allele_code_distribution(int32_t* const allele_codes, const int n_samples, const int n_alleles);
long test_write_unphased_pgen(
        const int32_t* const allele_codes,
        const int32_t allele_ct,
        const int pgen_file_mode,
        const long n_variants,
        const int n_samples,
        long &writtenVariantCount);
// integer constants to parallel PgenFileMode, for use when calling jni callable functions, which can't
// use the PgenFileMode enum provided by plink2
constexpr int PGEN_FILE_MODE_BACKWARD_SEEK = static_cast<int>(plink2::PgenWriteMode::kPgenWriteBackwardSeek);
constexpr int PGEN_FILE_MODE_WRITE_SEPARATE_INDEX = static_cast<int>(plink2::PgenWriteMode::kPgenWriteSeparateIndex);
constexpr int PGEN_FILE_MODE_WRITE_AND_COPY = static_cast<int>(plink2::PgenWriteMode::kPgenWriteAndCopy);

//******************* Tests *******************

// simple test to exercise throwing/catching of PgenException
BOOST_AUTO_TEST_CASE(test_exception_propagation) {
    const char* const expectedPropagationMessage = "Fake pgen exception";
    BOOST_REQUIRE_EXCEPTION(
            throw PgenException(expectedPropagationMessage),
            PgenException,
            [expectedPropagationMessage](PgenException ex) -> bool {
                return !strcmp(ex.what(), expectedPropagationMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(test_pglerr_string_conversion) {
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
constexpr boost::array<int, 3> s_pgenFileMode {
    PGEN_FILE_MODE_BACKWARD_SEEK,
    PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
    PGEN_FILE_MODE_WRITE_AND_COPY
};
BOOST_DATA_TEST_CASE(test_write_unphased_biallelic_small, s_pgenFileMode) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples
    constexpr int32_t allele_codes[] {0, 0, 0, 1, 1, 1 };
    long variantCount = 0L;
    // don't be fooled by the reference to the variable "sample" below; its a variable name introduced by
    // the boost macro to refer to the parameter for the test case, which in this test is the pgen file mode...
    const long file_size = test_write_unphased_pgen(allele_codes, 2,sample, n_variants, n_samples, variantCount);
    // only check for non-zero file size, since the file size varies with the file mode
    BOOST_REQUIRE_NE(file_size, 0);
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// verify the issue described here: https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ for bi-allelics
BOOST_AUTO_TEST_CASE(test_biallelic_one_allele_not_observed) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    constexpr int n_alleles = 2;
    // one variants's worth of allele codes - drawn from 2 alleles over 3 samples, but with only one
    // allele actually observed
    constexpr int32_t allele_codes[] {0, 0, 0, 0, 0, 0 };
    long variantCount = 0L;
    test_write_unphased_pgen(allele_codes, n_alleles, PGEN_FILE_MODE_WRITE_AND_COPY, n_variants, n_samples, variantCount);
}

// verify issue described here: https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ for multi-allelics
BOOST_AUTO_TEST_CASE(test_multiallelic_some_alleles_not_observed) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    constexpr int n_alleles = 7;
    // one variants's worth of allele codes - drawn from 7 alleles over 3 samples, but with only 3
    // alleles actually observed
    constexpr int32_t allele_codes[] {0, 1, 0, 5, 0, 4 };
    long variantCount = 0L;
    test_write_unphased_pgen(allele_codes, n_alleles, PGEN_FILE_MODE_WRITE_AND_COPY, n_variants, n_samples, variantCount);
}

// write a larger, bi-allelic pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(test_write_unphased_biallelic_large) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 2;
    // one variants's worth of allele codes drawn from 2 alleles
    int32_t *allele_codes = new int32_t[n_samples * 2];
    generate_allele_code_distribution(allele_codes, n_samples, n_alleles);
    long variantCount = 0L;
    const long file_size = test_write_unphased_pgen(allele_codes, n_alleles,PGEN_FILE_MODE_WRITE_AND_COPY, n_variants, n_samples, variantCount);
    delete[] allele_codes;

    BOOST_REQUIRE_EQUAL(file_size, 125450028); // cause thats what it is
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// write a larger, bi-allelic pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(test_write_unphased_multi_allelic_large) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 7;
    // synthesize one variants's worth of allele codes, with genotypes randomly drawn from 7 allele codes
    int32_t *allele_codes = new int32_t[n_samples * 2];
    generate_allele_code_distribution(allele_codes, n_samples, n_alleles);
    long variantCount = 0L;
    const long file_size = test_write_unphased_pgen(allele_codes, n_alleles,PGEN_FILE_MODE_WRITE_AND_COPY, n_variants, n_samples, variantCount);
    delete[] allele_codes;

    BOOST_REQUIRE_EQUAL(file_size, 911252530); // cause thats what it is
    BOOST_REQUIRE_EQUAL(variantCount, n_variants);
}

// claim that we're going to write 10 variants, but then don't write them
BOOST_AUTO_TEST_CASE(test_close_with_no_writes) {
    const char* const expectedNoWriteMessage = "closePgen called with number of variants written";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            10L,
            3,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            closePgen(pgenContext, 0),
            PgenException,
            [expectedNoWriteMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(test_invalid_allele_code) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples, with on bad allele code
    constexpr int32_t allele_codes[] {0, 0, 0, -17, 0, 0 };
    const char* const expectedMessageFragment = "Attempt to append invalid allele code";
    long variantCount;
    BOOST_REQUIRE_EXCEPTION(
            test_write_unphased_pgen(allele_codes, 2, PGEN_FILE_MODE_WRITE_AND_COPY, n_variants, n_samples, variantCount),
            PgenException,
            [expectedMessageFragment](PgenException ex) -> bool {
                return strstr(ex.what(), expectedMessageFragment);
            }
    );
}

// use a bogus pgen write mode
BOOST_AUTO_TEST_CASE(test_undefined_write_mode) {
    const char* const expectedInvalidModeMessage = "Invalid pgenWriteMode";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const int invalidWriteMode = -2; // must be one of 0, 1, 2
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, invalidWriteMode, 10L, 3, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

// use pgen write mode 0 (seek) with an unknown variant count
BOOST_AUTO_TEST_CASE(test_seek_write_mode_requires_known_variant_count) {
    const char* const expectedInvalidModeMessage = "requires a known variant count";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, PGEN_FILE_MODE_BACKWARD_SEEK, static_cast<long>(plink2::kPglMaxVariantCt), 3, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(test_invalid_sample_count) {
    const char* const expectedInvalidSampleCountMessage = "Invalid sample count";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const int invalidSampleCount = 0; // must be >= 1
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, PGEN_FILE_MODE_WRITE_AND_COPY, 10L, invalidSampleCount, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidSampleCountMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidSampleCountMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(test_invalid_variant_count) {
    const char* const expectedInvalidVariantCountMessage = "Invalid variant count";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const long invalidVariantCount = 0; // must be >= 1
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, PGEN_FILE_MODE_WRITE_AND_COPY, invalidVariantCount, 3, plink2::kPglMaxAltAlleleCt),
            PgenException,
            [expectedInvalidVariantCountMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidVariantCountMessage);
            }
    );
}

//******************* Local Test Utilities *******************

// generate allele codes for n_samples, distributed across values from n_alleles
void generate_allele_code_distribution(int32_t* const allele_codes, const int n_samples, const int n_alleles) {
    for (int i = 0; i < (n_samples * 2); i+=2) {
        allele_codes[i] = i % n_alleles;
        allele_codes[i+1] = (i + 1) % n_alleles;
    }
}

// the caller should call unlink() on the resulting file to cause it to be deleted
template<size_t N>
void createTempFile(const char* const nameTemplate, char (&outputFileName)[N]) {
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

// write a (unphased) pGEN file given allele codes (the same allele code vector is used for each variant),
// given # of variants, # of samples, and write mode
long test_write_unphased_pgen(
        const int32_t *allele_codes,
        const int32_t allele_ct,
        const int pgen_file_mode,
        const long n_variants,
        const int n_samples,
        long &variantCount) {
    char tmp_file_name[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmp_file_name);

    const pgenlib::PgenContext *const pgen_context = pgenlib::openPgen(
            tmp_file_name,
            pgen_file_mode,
            n_variants,
            n_samples,
            plink2::kPglMaxAltAlleleCt);
    BOOST_REQUIRE_NE(pgen_context, nullptr);

    for (int i = 0; i < n_variants; i++) {
        pgenlib::appendAlleles(pgen_context, allele_codes, allele_ct);
    }
    variantCount = getNumberOfVariantsWritten(pgen_context);
    closePgen(pgen_context, 0);

    // for now, just validate that the file has SOME contents; the enclosing pgen project has test code
    // that verifies the contents using plink2
    struct stat st;
    stat(tmp_file_name, &st);
    const long file_size = st.st_size;
    BOOST_REQUIRE_NE(file_size, 0);

    unlink(tmp_file_name);
    return file_size;
}
