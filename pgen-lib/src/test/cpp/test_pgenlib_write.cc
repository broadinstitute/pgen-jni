#define BOOST_TEST_MODULE pgen_write
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;

#include "pgenMeta.h"
#include "pgenIO.h"

BOOST_AUTO_TEST_CASE(test_open_pgen) {
    PGEN_META * const pgenMeta = openPgen("test_open_pgen.out", 10, 3);
    closePgen(pgenMeta);

    BOOST_TEST(true);
}
