#include "pgenContext.h"
#include "pgenException.h"
#include "pgenUtils.h"
#include "pgenlib_misc.h"

namespace pgenlib {

    PgenContext *openPgen(
            const char *cFilename,
            const long numberOfVariants,
            const long sampleCount) {
        PgenContext *const pGenContext = static_cast<PgenContext *const>(malloc(sizeof(PgenContext)));
        if (pGenContext == nullptr) {
            throw PgenException("Native code failure allocating PgenContext");
        }

        pGenContext->spgwp = static_cast<plink2::STPgenWriter *>(malloc(sizeof(plink2::STPgenWriter)));
        if (pGenContext->spgwp == nullptr) {
            free(pGenContext);
            throw PgenException("Native code failure allocating STPgenWriter");
        }

        uintptr_t alloc_cacheline_ct_ptr;
        uint32_t max_vrec_len;

        const plink2::PglErr init1Result = plink2::SpgwInitPhase1(cFilename, //filename
                                                                  nullptr,  // allele index offsets ( for multi allele)
                                                                  nullptr,  // non-ref flags
                                                                  (uint32_t) numberOfVariants, // number of variants
                                                                  (uint32_t) sampleCount, // sample count
                                                                  0, // optional max allele count
                                                                  plink2::kPgenWriteSeparateIndex, //todo PgenWriteMode == kPgenWriteSeparateIndex
                                                                  plink2::kfPgenGlobal0, //todo- is this right ? type: PgenGlobalFlags phase dosage gflags (genotype?)
                                                                  0, //  non-ref flags storage
                                                                  pGenContext->spgwp, // STPgenWriter * spgwp
                                                                  &alloc_cacheline_ct_ptr, //  uintptr_t* alloc_cacheline_ct_ptr
                                                                  &max_vrec_len);  // max vrec len ptr

        uint32_t bitvec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kBitsPerCacheline);

        //        if reterr != kPglRetSuccess:
        //            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
        //    throwIfPglErr(env, init1Result, "Initialization phase 1 failed");

        //        cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
        //        cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
        uint32_t genovec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kNypsPerCacheline);
        uint32_t dosage_main_cacheline_ct = plink2::DivUp(sampleCount, (2 * plink2::kInt32PerCacheline));

        unsigned char *spgw_alloc;
        if (plink2::cachealigned_malloc(
                (alloc_cacheline_ct_ptr + genovec_cacheline_ct + 3 * bitvec_cacheline_ct + dosage_main_cacheline_ct) *
                plink2::kCacheline, &spgw_alloc)) {
            throw PgenException("Native code failure allocating plink2::cachealigned_malloc");
        }

        //        SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
        //        self._genovec = <uintptr_t*>(&(spgw_alloc[alloc_cacheline_ct * kCacheline]))
        //        self._phasepresent = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct) * kCacheline]))
        //        self._phaseinfo = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + bitvec_cacheline_ct) * kCacheline]))
        //        self._dosage_present = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * kCacheline]))
        //        self._dosage_main = <uint16_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * kCacheline]))
        //        return
        SpgwInitPhase2(max_vrec_len, pGenContext->spgwp, spgw_alloc);
        pGenContext->genovec = (uintptr_t *) (&(spgw_alloc[alloc_cacheline_ct_ptr * plink2::kCacheline]));
        pGenContext->phasepresent = (uintptr_t *) (&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct) *
                                                                plink2::kCacheline]));
        pGenContext->phaseinfo = (uintptr_t *) (&(spgw_alloc[
                (alloc_cacheline_ct_ptr + genovec_cacheline_ct + bitvec_cacheline_ct) * plink2::kCacheline]));
        pGenContext->dosage_present = (uintptr_t *) (&(spgw_alloc[
                (alloc_cacheline_ct_ptr + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * plink2::kCacheline]));
        pGenContext->dosage_main = (uint16_t *) (&(spgw_alloc[
                (alloc_cacheline_ct_ptr + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * plink2::kCacheline]));

        return pGenContext;
    }

    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes ) {
        uintptr_t* genovec = pGenContext->genovec;
        plink2::AlleleCodesToGenoarrUnsafe(allele_codes, nullptr, plink2::SpgwGetSampleCt(pGenContext->spgwp), genovec, nullptr, nullptr);
        plink2::PglErr pglErr = plink2::SpgwAppendBiallelicGenovec(genovec, pGenContext->spgwp);
        throwOnPglErr(pglErr, "Native code failure adding genotypes");
    } 

    void closePgen(const PgenContext *const pGenContext) {
        const uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(pGenContext->spgwp);
        const uint32_t writtenVariantCt = plink2::SpgwGetVidx(pGenContext->spgwp);
        if (declaredVariantCt != writtenVariantCt) {
             snprintf(reservedForExceptionMessage,
                    kReservedMessageBufSize,
                    "PgenWriter.closePgen() called when number of written variants (%d) not equal to initially declared value (%d)",
                    writtenVariantCt,
                    declaredVariantCt);
            throw PgenException(reservedForExceptionMessage);
        }

        throwOnPglErr(SpgwFinish(pGenContext->spgwp), "Error in while closing pgen file");
    }

}