//

#ifndef PGEN_LIB_PGENCONTEXT_H
#define PGEN_LIB_PGENCONTEXT_H

#include "pgenlib_write.h"
#include "pgenlib_ffi_support.h"

namespace pgenlib {

    typedef struct PgenContext {
        plink2::STPgenWriter* spgwp;
        //uintptr_t alloc_cacheline_ct;
        uint32_t sampleCount;
        uint32_t max_vrec_len;
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

        //    cdef uintptr_t* _nonref_flags
        //    cdef PgenGlobalFlags _phase_dosage_gflags
        //    cdef uint32_t _allele_ct_limit

        //    # preallocate buffers we'll use repeatedly
        //    cdef uintptr_t* _genovec
        //    cdef uintptr_t* _patch_01_set
        //    cdef AlleleCode* _patch_01_vals
        //    cdef uintptr_t* _patch_10_set
        //    cdef AlleleCode* _patch_10_vals
        //    cdef uintptr_t* _phasepresent
        //    cdef uintptr_t* _phaseinfo
        //    cdef uintptr_t* _dosage_present
        //    cdef uint16_t* _dosage_main

    } PgenContext;
    
}
#endif //PGEN_LIB_PGENCONTEXT_H
