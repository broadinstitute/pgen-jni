//

#ifndef PGEN_LIB_PGENCONTEXT_H
#define PGEN_LIB_PGENCONTEXT_H

#include "pgenlib_write.h"
#include "pgenlib_ffi_support.h"

namespace pgenlib {

    typedef struct PgenContext {
        plink2::STPgenWriter* spgwp;
        // keep track of the arena memory so we can free it at the end
        unsigned char* spgw_alloc;
        // cdef uintptr_t* _nonref_flags
        // cdef PgenGlobalFlags _phase_dosage_gflags
        uint32_t allele_ct_limit;

        uintptr_t* genovec;         // genotype vector
        uintptr_t* patch_01_set;
        plink2::AlleleCode* patch_01_vals;
        uintptr_t* patch_10_set;
        plink2::AlleleCode* patch_10_vals;
        uintptr_t* phasepresent;
        uintptr_t* phaseinfo;
        uintptr_t* dosage_present;
        uint16_t* dosage_main;

        uint32_t sampleCount;
        uint32_t max_vrec_len;

        // non-plink2 fields
        uint32_t writeFlags; // keep track of whether the caller claims to have phasing data/multi-allelics
    } PgenContext;

}
#endif //PGEN_LIB_PGENCONTEXT_H
