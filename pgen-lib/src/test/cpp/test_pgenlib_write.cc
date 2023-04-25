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

//******************* Forward Declarations/Constants *******************
constexpr int TMP_FILENAME_SIZE = 4096;
template<size_t N> void createTempFile(const char* const nameTemplate, char (&outputFileName)[N]);
long test_write_pgen(
        const long n_variants,
        const int n_samples,
        const int pgen_file_mode,
        const int32_t* const allele_codes);
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

// write a small, bi-allelic pgen file once with each possible file write mode
static constexpr boost::array<int, 3> s_pgenFileMode {
    PGEN_FILE_MODE_BACKWARD_SEEK,
    PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
    PGEN_FILE_MODE_WRITE_AND_COPY
};
BOOST_DATA_TEST_CASE(test_write_biallelic_small, s_pgenFileMode) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples
    constexpr int32_t allele_codes[] {0, 0, 0, 0, 0, 0 };
    // don't be fooled by the reference to the variable "sample" below; its a variable name introduced by
    // the boost macro to refer to the parameter for the test case, which in this test is the pgen file mode...
    test_write_pgen(n_variants, n_samples, sample, allele_codes);
    // ignore the file size, since it varies with the file mode
}

// write a larger, bi-allelic pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(test_write_biallelic_large) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 2;
    // one variants's worth of allele codes - 2 alleles over n_samples
    int32_t *allele_codes = new int32_t[n_samples * 2];
    for (int i = 0; i < n_samples; i+=2) {
        allele_codes[i] = rand() % n_alleles;
        allele_codes[i+1] = rand() % n_alleles;
    }
    const long file_size = test_write_pgen(n_variants, n_samples, PGEN_FILE_MODE_WRITE_AND_COPY, allele_codes);
    delete[] allele_codes;

    //TODO: hm - for some reason, the file is 350028 on my Mac, but is 353340 on CI/linux
    //BOOST_REQUIRE_EQUAL(file_size, 350028); // cause thats what it is
}

// write a larger, bi-allelic pgen, using only file mode PGEN_FILE_MODE_WRITE_AND_COPY
BOOST_AUTO_TEST_CASE(test_write_multi_allelic_large) {
    constexpr long n_variants = 100000L;
    constexpr int n_samples = 10000;
    constexpr int n_alleles = 7;
    // one variants's worth of allele codes - 7 alleles over n_samples
    int32_t *allele_codes = new int32_t[n_samples * 2];
    for (int i = 0; i < n_samples; i+=2) {
        allele_codes[i] = rand() % n_alleles;
        allele_codes[i+1] = rand() % n_alleles;
    }
    const long file_size = test_write_pgen(n_variants, n_samples, PGEN_FILE_MODE_WRITE_AND_COPY, allele_codes);
    delete[] allele_codes;

    //TODO: hm - for some reason, the file is 350028 on my Mac, but is 353340 on CI/linux
    //BOOST_REQUIRE_EQUAL(file_size, 350028); // cause thats what it is
}

BOOST_AUTO_TEST_CASE(test_write_bad_allele_code) {
    constexpr long n_variants = 6;
    constexpr int n_samples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples, with on bad allele code
    constexpr int32_t allele_codes[] {0, 0, 0, -17, 0, 0 };
    //TODO: fix the error messages in the C++ code
    const char* const expectedMessage = "ConvertMultiAlleleCodesUnsafe";
    BOOST_REQUIRE_EXCEPTION(
            test_write_pgen(n_variants, n_samples, PGEN_FILE_MODE_WRITE_AND_COPY, allele_codes),
            PgenException,
            [expectedMessage](PgenException ex) -> bool {
                return strstr(ex.what(), expectedMessage);
            }
    );
}

// say we're going to write 10 variants, but don't write them
BOOST_AUTO_TEST_CASE(test_close_with_no_writes) {
    const char* const expectedNoWriteMessage = "number of written variants";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen(
            tmpFileName,
            PGEN_FILE_MODE_WRITE_AND_COPY,
            10L,
            3);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            closePgen(pgenContext),
            PgenException,
            [expectedNoWriteMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
}

// use a bogus pgen write mode
BOOST_AUTO_TEST_CASE(test_invalid_write_mode) {
    const char* const expectedInvalidModeMessage = "Invalid pgenWriteMode";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    const int invalidWriteMode = -2; // must be one of 0, 1, 2
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, invalidWriteMode, 10L, 3),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

//******************* Start Local Test Utils *******************

// writing of a pGEN file given allele codes (the same allele code vector is used for each variant),
// # variants, # samples, and write mode
long test_write_pgen(
        const long n_variants,
        const int n_samples,
        const int pgen_file_mode,
        const int32_t *allele_codes) {
    char tmp_file_name[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmp_file_name);

    const pgenlib::PgenContext *const pgen_context = pgenlib::openPgen(
            tmp_file_name,
            pgen_file_mode,
            n_variants,
            n_samples);
    BOOST_REQUIRE_NE(pgen_context, nullptr);
    for (int i = 0; i < n_variants; i++) {
        pgenlib::appendAlleles(pgen_context, allele_codes);
    }
    closePgen(pgen_context);

    // for now, just validate that the file has SOME contents; the enclosing pgen project has test code
    // that verifies the contents using plink2
    struct stat st;
    stat(tmp_file_name, &st);
    const long file_size = st.st_size;
    BOOST_REQUIRE_NE(file_size, 0);

    unlink(tmp_file_name);
    return file_size;
}

// the caller should call unlink() on the resulting file to cause it to be deleted (careful - if its open
// it will be deleted on close)
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
    // the file is open, so close it, but DON'T call unlink for it to be deleted - leave that to the caller...
    close(fDesc);
}
