#define BOOST_TEST_MODULE pgen_write
#include <boost/test/included/unit_test.hpp>
#include "pgenException.h"

#include "pgenContext.h"
#include "pgenIO.h"

using namespace boost::unit_test;
using namespace pgenlib;

const char *expectedPropagationMessage = "Fake pgen exception";
bool validatePropagationMessage(PgenException ex) {
    return !strcmp(ex.what(), expectedPropagationMessage);
}
BOOST_AUTO_TEST_CASE(test_pgen_exception_propagation) {
    BOOST_REQUIRE_EXCEPTION(
            throw PgenException(expectedPropagationMessage),
            PgenException,
            ::validatePropagationMessage
    );
};

//BOOST_AUTO_TEST_CASE(test_write_pgen) {
//    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen("test_open.pgen", 10, 3);
//    BOOST_REQUIRE_NE(pgenContext, nullptr);
//    //TODO: appendAlleles
//    closePgen(pgenContext);
//    //TODO: check for the file
//}

const char *expectedNoWriteMessage = "number of written variants";
bool validateNoWriteMessage(PgenException ex) {
    return strstr(ex.what(), expectedNoWriteMessage);
}
BOOST_AUTO_TEST_CASE(test_close_pgen_with_no_writes) {
    // say we're going to write 10 variants, but don't write them
    const pgenlib::PgenContext *const pgenContext = pgenlib::openPgen("test_open.pgen", 10, 3);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    BOOST_REQUIRE_EXCEPTION(
            closePgen(pgenContext),
            PgenException,
            ::validateNoWriteMessage
    );
}
