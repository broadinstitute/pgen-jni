#define BOOST_TEST_MODULE pgen_write
#include <sys/stat.h>
#include <stdio.h>

#include <boost/test/included/unit_test.hpp>
#include "pgenException.h"
#include "pgenContext.h"
#include "pgenIO.h"

using namespace boost::unit_test;
using namespace pgenlib;

const int TMP_FILENAME_SIZE = 4096;
template<size_t N> void createTempFile(const char *nameTemplate, char (&oututFileName)[N]);

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

BOOST_AUTO_TEST_CASE(test_write_pgen) {
    long numberOfVariants = 6L;
    int numberOfSamples = 3;
    // one variants's worth of allele codes - 2 alleles over 3 samples
    int32_t allele_codes[] {0, 0, 0, 0, 0, 0 };

    char tmpFileName[TMP_FILENAME_SIZE];
    createTempFile("test_write.pgen", tmpFileName);

    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen(
            tmpFileName,
            2,
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
    long fileSize = st.st_size;
    BOOST_REQUIRE_NE(fileSize, 0);

    unlink(tmpFileName);
}

// say we're going to write 10 variants, but don't write them
BOOST_AUTO_TEST_CASE(test_close_pgen_with_no_writes) {
    const char *expectedNoWriteMessage = "number of written variants";
    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen("test_open.pgen", 1, 10L, 3);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    BOOST_REQUIRE_EXCEPTION(
            // note, this will throw and no cleanup will have been done...
            closePgen(pgenContext),
            PgenException,
            [expectedNoWriteMessage](PgenException ex) -> bool  {
                return strstr(ex.what(), expectedNoWriteMessage);
            }
    );
}

// call unlink(nameBuff) on thr resulting file to cause (careful - if its open it will be deleted on close)
template<size_t N>
void createTempFile(
        const char *nameTemplate,
        char (&oututFileName)[N]) // caller allocates and frees
{
    //this is deprecated (and maybe a little sketchy), but works nicely to obtain a tmp dir location
    std::string tmpPath = std::tmpnam(nullptr);
    snprintf(oututFileName, N, "%s_pgenBoostXXXXXX%s", tmpPath.c_str(), nameTemplate);

    // we don't actually need to create the file here, just reserve it
    int fDesc = mkstemps(oututFileName, strlen(nameTemplate));
    if (fDesc < 1) {
        char errMessage[pgenlib::kReservedMessageBufSize];
        snprintf(errMessage,
                 pgenlib::kReservedMessageBufSize,
                 "Temp file creation failed for (%s) with error(%s)",
                 oututFileName,
                 strerror(errno));
        throw PgenException(errMessage);
    }
    // the file is open, so close it, but DON'T call unlink until we're ready for it to be deleted...
    close(fDesc);
}
