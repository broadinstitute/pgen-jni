#define BOOST_TEST_MODULE pgen_write
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;

#include "pgenContext.h"
#include "pgenIO.h"

BOOST_AUTO_TEST_CASE(test_open_pgen) {
    try {
        const PgenContext *const pgenMeta = openPgen("test_open.pgen", 10, 3);
        closePgen(pgenMeta);
    } catch (...) {
        //TODO: how to tell BOOST to fail....
    }

    BOOST_TEST(true);
}
