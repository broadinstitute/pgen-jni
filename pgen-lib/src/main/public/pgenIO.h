//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

namespace pgenlib {

    // sentinel to use as a variant count to signal that the variant count is unknown
    const uint32_t kVariantCountUnknown = plink2::kPglMaxVariantCt;

    // the public interface to the pgen writer
    PgenContext *openPgen(const char *cFilename, const int pgenWriteModeInt, const unsigned int pgenWriteFlags, const long variantCount, const int sampleCount, const int maxAltAlleles);
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes, const unsigned char* phase_bytes, int32_t allele_ct);
    void closePgen(const PgenContext *const pGenContext, const long nDroppedVariants);
    long getNumberOfVariantsWritten(const PgenContext *const pGenContext);

}
#endif //PGEN_LIB_PGENIO_H
