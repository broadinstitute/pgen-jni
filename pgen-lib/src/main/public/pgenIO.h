//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

namespace pgenlib {

    PgenContext *openPgen(const char *cFilename, const long numberOfVariants, const long sampleCount);
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes );
    void closePgen(const PgenContext *const pGenContext);
}
#endif //PGEN_LIB_PGENIO_H
