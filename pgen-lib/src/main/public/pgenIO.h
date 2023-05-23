//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

namespace pgenlib {

    // calleable via jni
    PgenContext *openPgen(const char *cFilename, const int pgenWriteModeInt, const long variantCount, const int sampleCount, const int maxAltAlleles);
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes, int32_t allele_ct);
    void closePgen(const PgenContext *const pGenContext, const long nDroppedVariants);
    long getNumberOfVariantsWritten(const PgenContext *const pGenContext);

}
#endif //PGEN_LIB_PGENIO_H
