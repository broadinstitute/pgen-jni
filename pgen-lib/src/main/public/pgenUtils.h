//
#ifndef PGEN_LIB_PGENUTILS_H
#define PGEN_LIB_PGENUTILS_H
#include "pgenException.h"
#include "plink2_base.h"

namespace pgenlib {

    // Throw an exception for an error that originated in underlying *PLINK* code.
    void throwOnPglErr(const plink2::PglErr result, const char* message) {
        if (result) { // PglErr evaluates to true if there was an error
            //TODO: translate the PglErr enum value to a string
            const int MSG_BUF_SIZE = 4096;
            char msg[MSG_BUF_SIZE];
            snprintf(msg, MSG_BUF_SIZE, "%s (PglErr: %d)", message, int32_t(result));
            throw PgenException(msg);
        }
    }
 }

#endif //PGEN_LIB_PGENUTILS_H
