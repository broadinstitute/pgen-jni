#define BOOST_TEST_MODULE pgen_write
#include <boost/test/included/unit_test.hpp>
#include "pgenException.h"

#include "pgenContext.h"
#include "pgenIO.h"

using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(test_open_pgen) {
    const PgenContext *const pgenContext = openPgen("test_open.pgen", 10, 3);
    BOOST_REQUIRE_NE(pgenContext, nullptr);
    closePgen(pgenContext);
}

const char *expectedMessage = "Fake pgen exception";
bool validateExceptionMessage(PgenException ex) {
    return !strcmp(ex.what(), expectedMessage);
}
BOOST_AUTO_TEST_CASE(test_pgen_exception_propagation) {
    BOOST_REQUIRE_EXCEPTION(
            throw PgenException(expectedMessage),
            PgenException,
            ::validateExceptionMessage
    );
};
