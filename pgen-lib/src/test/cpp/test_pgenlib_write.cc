#define BOOST_TEST_MODULE pgen_write
#include <boost/test/included/unit_test.hpp>
#include "pgenException.h"

#include "pgenContext.h"
#include "pgenIO.h"

using namespace boost::unit_test;

//TODO: how to tell BOOST the test failed....
//TODO: BOOST - failure on exception ?

BOOST_AUTO_TEST_CASE(test_open_pgen) {
    const PgenContext *const pgenContext = openPgen("test_open.pgen", 10, 3);
    BOOST_TEST(pgenContext != nullptr);
    closePgen(pgenContext);
}

BOOST_AUTO_TEST_CASE(test_pgen_exception_propagation) {
    char *expectedMessage = "Fake pgen exception";
    bool isExpectedMessage = true;

    try {
        throw PgenException(expectedMessage);
    } catch (PgenException& e) {
        isExpectedMessage = !strcmp(e.what(), expectedMessage);
    }
    BOOST_TEST(isExpectedMessage);
}
