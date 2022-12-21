//

#ifndef PGEN_LIB_PGENMETA_H
#define PGEN_LIB_PGENMETA_H

#include "pgenlib_write.h"
#include "pgenlib_ffi_support.h"

typedef struct PGEN_META {
    plink2::STPgenWriter* spgwp;
    uintptr_t alloc_cacheline_ct_ptr;
    uint32_t max_vrec_len;

    uintptr_t* genovec; //Genotype vector
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    uintptr_t* dosage_present;
    uint16_t* dosage_main;

} pgenMeta;

#endif //PGEN_LIB_PGENMETA_H
