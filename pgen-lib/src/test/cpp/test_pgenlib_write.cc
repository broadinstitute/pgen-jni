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

const int TMP_FILENAME_SIZE = 4096;
template<size_t N> void createTempFile(const char *nameTemplate, char (&outputFileName)[N]);

// simple test to exercise throwing/catching of PgenException
BOOST_AUTO_TEST_CASE(test_pgen_exception_propagation) {
    const char *expectedPropagationMessage = "Fake pgen exception";
    BOOST_REQUIRE_EXCEPTION(
            throw PgenException(expectedPropagationMessage),
            PgenException,
            [expectedPropagationMessage](PgenException ex) -> bool {
                return !strcmp(ex.what(), expectedPropagationMessage);
            }
    );
}

BOOST_AUTO_TEST_CASE(test_pgl_string_conversion) {
    const char *expectedMessage = "kPglRetNotYetSupported";
    BOOST_REQUIRE_EXCEPTION(
            throwOnPglErr(plink2::PglErr::ec::kPglRetNotYetSupported, "Testing PglErr conversion"),
            PgenException,
            [expectedMessage](PgenException ex) -> bool {
                return strstr(ex.what(), expectedMessage);
            }
    );
}

// simple writing of a pGEN file with each possible file write mode
static const boost::array<int, 3> s_pgenFileMode { 0, 1, 2 };
BOOST_DATA_TEST_CASE(test_write_pgen, s_pgenFileMode) {
    const long numberOfVariants = 6L;
    const int numberOfSamples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples
    const int32_t allele_codes[] {0, 0, 0, 0, 0, 0 };

    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);

    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen(
            tmpFileName,
            sample, // the variable "sample" is magically introduced by the BOOST_DATA_TEST_CASE macro
            numberOfVariants,
            numberOfSamples);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    for (int i = 0; i < numberOfVariants; i++) {
        pgenlib::appendAlleles(pgenContext, allele_codes);
    }
    closePgen(pgenContext);

    // for now, just validate that the file has SOME contents; the enclosing pgen project has test code
    // that verifies the contents using plink2
    struct stat st;
    stat(tmpFileName, &st);
    const long fileSize = st.st_size;
    BOOST_REQUIRE_NE(fileSize, 0);

    unlink(tmpFileName);
}

// say we're going to write 10 variants, but don't write them
BOOST_AUTO_TEST_CASE(test_close_pgen_with_no_writes) {
    const char *expectedNoWriteMessage = "number of written variants";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen(tmpFileName, 1, 10L, 3);
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

BOOST_AUTO_TEST_CASE(test_invalid_write_mode) {
    const char *expectedInvalidModeMessage = "Invalid pgenWriteMode";
    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);
    unlink(tmpFileName);
    BOOST_REQUIRE_EXCEPTION(
            pgenlib::openPgen(tmpFileName, -2, 10L, 3),
            PgenException,
            [expectedInvalidModeMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedInvalidModeMessage);
            }
    );
}

// the caller should call unlink() on the resulting file to cause it to be deleted (careful - if its open
// it will be deleted on close)
template<size_t N>
void createTempFile(const char *nameTemplate, char (&outputFileName)[N]) {
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
