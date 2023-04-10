//
#ifndef PGEN_LIB_PGENUTILS_H
#define PGEN_LIB_PGENUTILS_H

#include "pgenException.h"
#include "plink2_base.h"

namespace pgenlib {
    // Throw an exception for an error that originated in underlying *PLINK* code.
    void throwOnPglErr(const plink2::PglErr result, const char* message);
 }

#endif //PGEN_LIB_PGENUTILS_H
