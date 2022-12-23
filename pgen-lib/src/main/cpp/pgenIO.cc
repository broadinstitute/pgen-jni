#include "pgenMeta.h"
#include "pgenlib_misc.h"

PGEN_META *openPgen (
        const char* cFilename,
        long numberOfVariants,
        long sampleCount) {
    PGEN_META * const pGenMeta = (PGEN_META*) malloc(sizeof(PGEN_META));
    if (pGenMeta == NULL){
        //todo
    }

    pGenMeta->spgwp = (plink2::STPgenWriter*) malloc(sizeof(plink2::STPgenWriter));
    if (pGenMeta->spgwp == NULL){
        //todo
    }

    uintptr_t alloc_cacheline_ct_ptr;
    uint32_t max_vrec_len;

    const plink2::PglErr init1Result = plink2::SpgwInitPhase1(  cFilename, //filename
                                                                NULL,  // allele index offsets ( for multi allele)
                                                                NULL,  // non-ref flags
                                                                (uint32_t) numberOfVariants, // number of variants
                                                                (uint32_t) sampleCount, // sample count
                                                                0, // optional max allele count
                                                                plink2::kPgenWriteSeparateIndex, //todo PgenWriteMode == kPgenWriteSeparateIndex
                                                                plink2::kfPgenGlobal0, //todo- is this right ? type: PgenGlobalFlags phase dosage gflags (genotype?)
                                                                0, //  non-ref flags storage
                                                                pGenMeta->spgwp , // STPgenWriter * spgwp
                                                                &alloc_cacheline_ct_ptr , //  uintptr_t* alloc_cacheline_ct_ptr
                                                                &max_vrec_len);  // max vrec len ptr

    uint32_t bitvec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kBitsPerCacheline);

    //        if reterr != kPglRetSuccess:
    //            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//    checkPglErr(env, init1Result, "Initialization phase 1 failed");

    //        cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
    //        cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
    uint32_t genovec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kNypsPerCacheline);
    uint32_t dosage_main_cacheline_ct = plink2::DivUp(sampleCount, (2 * plink2::kInt32PerCacheline));

    unsigned char* spgw_alloc;
    if (plink2::cachealigned_malloc((alloc_cacheline_ct_ptr + genovec_cacheline_ct + 3 * bitvec_cacheline_ct + dosage_main_cacheline_ct) * plink2::kCacheline, &spgw_alloc)){
        //todo handle malloc fail
    }

    //        SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
    //        self._genovec = <uintptr_t*>(&(spgw_alloc[alloc_cacheline_ct * kCacheline]))
    //        self._phasepresent = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct) * kCacheline]))
    //        self._phaseinfo = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + bitvec_cacheline_ct) * kCacheline]))
    //        self._dosage_present = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * kCacheline]))
    //        self._dosage_main = <uint16_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * kCacheline]))
    //        return
    SpgwInitPhase2(max_vrec_len, pGenMeta->spgwp, spgw_alloc);
    pGenMeta->genovec = (uintptr_t*)(&(spgw_alloc[alloc_cacheline_ct_ptr * plink2::kCacheline]));
    pGenMeta->phasepresent = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct) * plink2::kCacheline]));
    pGenMeta->phaseinfo = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + bitvec_cacheline_ct) * plink2::kCacheline]));
    pGenMeta->dosage_present = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * plink2::kCacheline]));
    pGenMeta->dosage_main = (uint16_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * plink2::kCacheline]));

    return pGenMeta;
}

void closePgen(PGEN_META * const pGenMeta) {
    const uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(pGenMeta->spgwp);
    const uint32_t writtenVariantCt = plink2::SpgwGetVidx(pGenMeta->spgwp);
    if ( declaredVariantCt != writtenVariantCt) {
//        char errbuff[200];
//        sprintf(errbuff, "PgenWriter.close() called when number of written variants (%d) unequal to initially declared value (%d)", writtenVariantCt, declaredVariantCt);
//        jclass exceptionClass = env->FindClass("org/broadinstitute/pgen/PgenJniException");
//        env->ThrowNew(exceptionClass, errbuff);
        return;
    }

    //checkPglErr(env, SpgwFinish(bookKeepingPtr->spgwp), "Error while closing PgenFile");
}
