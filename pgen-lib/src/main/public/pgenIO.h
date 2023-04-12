//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

namespace pgenlib {

    // calleable via jni
    PgenContext *openPgen(const char *cFilename, const int pgenWriteModeInt, const long numberOfVariants, const int sampleCount);
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes );
    void closePgen(const PgenContext *const pGenContext);
}
#endif //PGEN_LIB_PGENIO_H
