//

#ifndef PGEN_LIB_PGENIO_H
#define PGEN_LIB_PGENIO_H
#include "pgenContext.h"

// the public interface to the PGEN writer
namespace pgenlib {

    // sentinel to use as a variant count to signal that the variant count is unknown
    constexpr uint32_t kVariantCountUnknown = plink2::kPglMaxVariantCt;

    // write flag values
    constexpr uint32_t kWriteFlagPreservePhasing = 0x1;
    constexpr uint32_t kWriteFlagMultiAllelic = 0x2;

    PgenContext *openPgen(const char *cFilename, const uint32_t pgenWriteModeInt, const uint32_t pgenWriteFlags, const long variantCount, const int sampleCount, const int maxAltAlleles);
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes, const unsigned char* phase_bytes, int32_t allele_ct);
    long getNumberOfVariantsWritten(const PgenContext *const pGenContext);
    void closePgen(const PgenContext *const pGenContext, const long nDroppedVariants);

}
#endif //PGEN_LIB_PGENIO_H
