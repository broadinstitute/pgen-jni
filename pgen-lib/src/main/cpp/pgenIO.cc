#include <cmath>

using namespace std;

#include "pgenContext.h"
#include "pgenException.h"
#include "pgenUtils.h"
#include "pgenlib_misc.h"
#include "pgenlib_write.h"

namespace pgenlib {

    plink2::PgenWriteMode validatePgenWriteMode(const uint32_t anInt);
    PgenContext *initPgenContext(
            const char* cFilename,
            const plink2::PgenWriteMode pgenWriteMode,
            const long variantCount,
            const int sampleCount,
            const int maxAltAlleles);
    static const int kErrMessageBufSize = 1024;

    /**
     * Start a new PGEN write session, and return a pointer to a PgenContext for the writer.
     *
     * The PgenContext can be used to write an entire pgen file. depending on
     * the pgen file mode used (see below) a .pgen.pgi file may also be created.
     *
     * Only diploid genomes are supported.
     *
     * The allele codes for the genotypes for each variant can be provided by passing the PgenContext returned by
     * this function to a series of calls to appendAlleles, after which the PgenContext session should be closed
     * via a call to closePgen.
     *
     * An example pgen writer lifecyle is illustrated here:
     *
     *      const pgenlib::PgenContext *const pgen_context = pgenlib::openPgen(
     *          file_name,
     *          pgen_file_mode,
     *          n_variants,
     *          n_samples,
     *          plink2::kPglMaxAltAlleleCt);
     *
     *      for (int i = 0; i < n_variants; i++) {
     *          pgenlib::appendAlleles(pgen_context, allele_codes, allele_ct);
     *      }
     *      long variantCount = getNumberOfVariantsWritten(pgen_context);
     *      closePgen(pgen_context, 0);
     *
     *  Once the PgenContext has been closed, it can no longer be used to write allele codes.
     *
     * @param cFilename - the pgen file to write
     * @param pgenWriteModeInt - an integer representing the file mode, with permitted values drawn from integer
     * values of plink2::PgenWriteMode (1, 2 or 3). An exception will be thrown if any other value is provided. this
     * determines the pgen file mode that is used (i.e, whether there is a separate .pgi index)
     * @param variantCount - the number of variant to be written. if fewer variants are written, an exception will
     * be thrown when the writer is closed by a call to closePgen. must be in the range 1..plink2::kPglMaxVariantCt
     * @param sampleCount - the number of samples (genotypes) in the data set. Must be > 0.
     * @param maxAltAlleles - the maximum number of alleles for any variant that will be written - this determines the
     * range of valid (zero based) allele codes that can be provided when writing genotypes to this writer. Must be
     * in the range 2..plink2::kPglMaxAltAlleleCt
     *
     * @return a PgenContext
     */
    PgenContext *openPgen(
            const char* cFilename,
            const int pgenWriteModeInt,
            const long variantCount,
            const int sampleCount,
            const int maxAltAlleles) {

        // validate the requested pgen write mode, and sample and variant counts
        plink2::PgenWriteMode pgenWriteMode = validatePgenWriteMode(pgenWriteModeInt);
        if (sampleCount < 1) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid sample count: %d. At least 1 sample is required.",
                     sampleCount);
            throw PgenException(errMessageBuff); // PgenException makes a copy of errMessageBuff
        } else if (variantCount < 1) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid variant count: %ld. Variant count must be > 0.",
                     variantCount);
            throw PgenException(errMessageBuff);  // PgenException makes a copy of errMessageBuff
        } else if (variantCount > plink2::kPglMaxVariantCt) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid variant count: %ld exceeds maximum allowable variant count: %d.",
                     variantCount,
                     plink2::kPglMaxVariantCt);
            throw PgenException(errMessageBuff); // PgenException makes a copy of errMessageBuff
        } else if (maxAltAlleles < 2) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid max alt allele count: %d must be at least 2.",
                     maxAltAlleles);
            throw PgenException(errMessageBuff);  // PgenException makes a copy of errMessageBuff
        } else if (maxAltAlleles > plink2::kPglMaxAltAlleleCt) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid max alt allele count: %d exceeds maximum allowable alt allele count: %d.",
                     maxAltAlleles,
                     plink2::kPglMaxAltAlleleCt);
            throw PgenException(errMessageBuff);  // PgenException makes a copy of errMessageBuff
        }

        return initPgenContext(cFilename, pgenWriteMode, variantCount, sampleCount, maxAltAlleles);
    }

    PgenContext *initPgenContext(
            const char* cFilename,
            const plink2::PgenWriteMode pgenWriteMode,
            const long variantCount,
            const int sampleCount,
            const int maxAltAlleles) {

        PgenContext* pGenContext = static_cast<PgenContext *const>(malloc(sizeof(PgenContext)));
        if (pGenContext == nullptr) {
            throw PgenException("Native code failure allocating PgenContext");
        }
        pGenContext->spgwp = static_cast<plink2::STPgenWriter *>(malloc(sizeof(plink2::STPgenWriter)));
        if (pGenContext->spgwp == nullptr) {
            free(pGenContext);
            throw PgenException("Native code failure allocating STPgenWriter");
        }

        // convert sampleCount and variantCount to the types plink uses
        pGenContext->sampleCount = static_cast<uint32_t>(sampleCount);
        uint32_t variant_ct = static_cast<uint32_t>(variantCount);
        // total max allele count is max alt count + 1
        pGenContext->allele_ct_limit = static_cast<uint32_t>(maxAltAlleles + 1);

        //TODO: what is non_ref_strage_flags values 1, 2, 3 ?

        uint32_t bitvec_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, plink2::kBitsPerCacheline);
        uintptr_t alloc_cacheline_ct = 0;
        const plink2::PglErr init1Result = plink2::SpgwInitPhase1(cFilename, //filename
                                                                nullptr,  // allele index offsets ( for multi allele)
                                                                nullptr,  // non-ref flags
                                                                variant_ct, // number of variants
                                                                pGenContext->sampleCount, // sample count
                                                                pGenContext->allele_ct_limit, // optional max allele count
                                                                pgenWriteMode,
                                                                plink2::kfPgenGlobal0,
                                                                1, //  non-ref flags storage
                                                                pGenContext->spgwp, // STPgenWriter * spgwp
                                                                &alloc_cacheline_ct, //  uintptr_t* alloc_cacheline_ct
                                                                &pGenContext->max_vrec_len);
        throwOnPglErr(init1Result, "plink2 initialization (SpgwInitPhase1 failed)");

        uint32_t genovec_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, plink2::kNypsPerCacheline);
        uint32_t patch_01_vals_cacheline_ct = plink2::DivUp(pGenContext->sampleCount * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t patch_10_vals_cacheline_ct = plink2::DivUp(pGenContext->sampleCount * 2 * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t dosage_main_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, (2 * plink2::kInt32PerCacheline));

        // There are two copies of pgenlib.pyx in the plink2 build, that have many differences. One uses + 3 for
        // this calculation, and one uses + 5. Prefer the one in src, and go with + 5.
        // keep the pointer to the arena block in pGenContext so we can free it at the end
        if (plink2::cachealigned_malloc(
                (alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) *
                plink2::kCacheline, &pGenContext->spgw_alloc)) {
            throw PgenException("Native code failure (cachealigned_malloc) allocating spgw_alloc");
        }
        SpgwInitPhase2(pGenContext->max_vrec_len, pGenContext->spgwp, pGenContext->spgw_alloc);

        unsigned char* spgw_alloc_iter = &(pGenContext->spgw_alloc[alloc_cacheline_ct * plink2::kCacheline]);

        pGenContext->genovec = (uintptr_t *) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_01_set = (uintptr_t*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_01_vals = (plink2::AlleleCode*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_10_set = (uintptr_t*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_10_vals = (plink2::AlleleCode*) spgw_alloc_iter;

        // we're not using the phase info yet
        spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * plink2::kCacheline]);
        pGenContext->phasepresent = (uintptr_t *) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->phaseinfo = (uintptr_t*) spgw_alloc_iter;

        // we probably don't need these dosage blocks, but per this comment from the python
        // code it probably doesnt much matter:
        // # Could skimp on dosage/phase, but that doesn't gain us much.
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->dosage_present = (uintptr_t*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->dosage_main = (uint16_t*) spgw_alloc_iter;

        //    # bugfix (16 Apr 2023): SpgwAppendBiallelicGenovec[Hphase] assumes
        //    # trailing bits are clear
        //    self._genovec[(sample_ct - 1) // kBitsPerWordD2] = 0
        pGenContext->genovec[(pGenContext->sampleCount - 1) / plink2::kBitsPerWordD2] = 0; // floor division
        pGenContext->phasepresent[(pGenContext->sampleCount - 1) / plink2::kBitsPerWord] = 0; // floor division

        return pGenContext;
    }

    /**
     * Append on variant's worth of allele code (genotypes) to a pgen file.
     *
     * @param pGenContext - the PgenContext for the writer
     * @param allele_codes - array of allele codes to be written
     * @param allele_ct - the number of possible allele values for this variant (not the number of unique alleles
     * that are ACTUALLY observed/present in allele_codes)
     */
    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes, const int32_t allele_ct) {
        //TODO: for now just declare allPhased == false since there is no phaseinfo provided;
        bool allPhased = false;
        uintptr_t *phaseinfo = nullptr;
        uint32_t patch_01_ct;
        uint32_t patch_10_ct;
        int32_t observed_allele_ct = plink2::ConvertMultiAlleleCodesUnsafe(
                allele_codes,
                nullptr,
                pGenContext->sampleCount,
                pGenContext->genovec,
                pGenContext->patch_01_set,
                pGenContext->patch_01_vals,
                pGenContext->patch_10_set,
                pGenContext->patch_10_vals,
                &patch_01_ct,
                &patch_10_ct,
                nullptr,
                phaseinfo);
        if (observed_allele_ct == -1) {
            // it would be nice if we could determine what the invalid code is
            throw PgenException("Attempt to append invalid allele code (plink2::ConvertMultiAlleleCodesUnsafe)");
        }
        uint32_t write_allele_ct = static_cast<uint32_t>(observed_allele_ct);
        uint32_t unsigned_allele_ct = static_cast<uint32_t>(allele_ct);
        if (write_allele_ct > pGenContext->allele_ct_limit) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "plink2::ConvertMultiAlleleCodesUnsafe found more allele codes (%u) than specified in allele_ct_limit (%u); you may need to construct the PgenWriter with a higher allele_ct_limit setting",
                     write_allele_ct,
                     pGenContext->allele_ct_limit);
            throw PgenException(errMessageBuff);
        }
        if (unsigned_allele_ct < write_allele_ct) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "plink2::ConvertMultiAlleleCodesUnsafe called with more alleles (%u) than stated in allele_ct (%u)",
                     write_allele_ct,
                     unsigned_allele_ct);
            throw PgenException(errMessageBuff);
        } else if (unsigned_allele_ct > pGenContext->allele_ct_limit) {
            // hm, this branch is actually not dependent on the call to ConvertMultiAlleleCodesUnsafe, and could
            // be done right at the start of the function, but I'll keep it here to match the flow of the python code
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "plink2::ConvertMultiAlleleCodesUnsafe called with allele_ct (%u) > allele_ct_limit (%u)",
                     unsigned_allele_ct,
                     pGenContext->allele_ct_limit);
            throw PgenException(errMessageBuff);
        }
        write_allele_ct = unsigned_allele_ct;

        //cout << " writing allele_ct " << allele_ct << " obs ct " << observed_allele_ct << "\n";
        //cout.flush();
        
        plink2::PglErr pglErr;
        if (!allPhased) {
            if ((patch_01_ct == 0) and (patch_10_ct == 0)) {
                // if (allele_ct != 2) {
                //     cout << "Writing multi-allelic " << allele_ct << " as bi-allelic";
                // }
                pglErr = SpgwAppendBiallelicGenovec(pGenContext->genovec, pGenContext->spgwp);
            } else {
                // if (allele_ct == 2) {
                //     cout << "Writing bi-allelic " << allele_ct << " as multi-allelic";
                // }
                pglErr = SpgwAppendMultiallelicSparse(
                        pGenContext->genovec,
                        pGenContext->patch_01_set,
                        pGenContext->patch_01_vals,
                        pGenContext->patch_10_set,
                        pGenContext->patch_10_vals,
                        write_allele_ct,
                        patch_01_ct,
                        patch_10_ct,
                        pGenContext->spgwp);
            }
        } else {
          if ((patch_01_ct == 0) and (patch_10_ct == 0)) {
              pglErr = SpgwAppendBiallelicGenovecHphase(
                      pGenContext->genovec,
                      NULL,
                      phaseinfo,
                      pGenContext->spgwp);
          } else {
              pglErr = SpgwAppendMultiallelicGenovecHphase(
                      pGenContext->genovec,
                      pGenContext->patch_01_set,
                      pGenContext->patch_01_vals,
                      pGenContext->patch_10_set,
                      pGenContext->patch_10_vals,
                      NULL,
                      phaseinfo,
                      write_allele_ct,
                      patch_01_ct,
                      patch_10_ct,
                      pGenContext->spgwp);
          }
        }
        throwOnPglErr(pglErr, "appendAlleles");
    }

    /**
     * Close a pGenContext, flush the output, and close the pgen file.
     *
     * @param pGenContext - the pgen context for this writer
     * @param numVariantsDropped - the number of variants dropped (the number f variants dropped, plus the number
     * of variants written, must equal the number of variants projected to be written when the pgen context was
     * initially opened. otherwise a pgenlib::PGenException will be thrown.
     */
    void closePgen(const PgenContext *const pGenContext, const long numVariantsDropped) {
        const uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(pGenContext->spgwp);
        const uint32_t writtenVariantCt = plink2::SpgwGetVidx(pGenContext->spgwp);
        const uint32_t droppedVariantCt = static_cast<uint32_t>(numVariantsDropped);

        if ((declaredVariantCt - droppedVariantCt) != writtenVariantCt) {
            //TODO: the plink2 python implementation throws on close if you haven't written as many variants as you
            // initially claimed you would, so we do too (after accounting for variants dropped due to exceeding the
            // maximum allele threshold). But, we've come this far - do we REALLY want to throw now ???
            char errMessage[kErrMessageBufSize];
            snprintf(reservedForExceptionMessage,
                     kErrMessageBufSize,
                      "closePgen called with number of variants written (%d) not equal to (declared - dropped)  (%d - %d = %d)",
                      writtenVariantCt,
                      declaredVariantCt,
                      droppedVariantCt,
                      declaredVariantCt - droppedVariantCt);
            throw PgenException(reservedForExceptionMessage);
        }

        throwOnPglErr(SpgwFinish(pGenContext->spgwp), "Error closing pgen file: SpgwFinish");

        //TODO: there might be a bug in plink pgen-lib, since I think the plink2 VCF importer only does one or the
        // other of SpgwFinish and CleanupSpgw (SpgwFinish on success, CleanupSpgw on failure), but not both. but
        // if we don't do both, the output doesn't seem to get flushed...
        plink2::PglErr cleanupErr;
        plink2::BoolErr bErr = CleanupSpgw(pGenContext->spgwp, &cleanupErr);
        if (bErr) {
            throwOnPglErr(cleanupErr, "Error cleaning up on pgen close: CleanupSpgw");
        }
        free(pGenContext->spgwp);
        plink2::aligned_free(pGenContext->spgw_alloc);
        free(reinterpret_cast<void*>(const_cast<PgenContext *>(pGenContext)));
    }

    long getNumberOfVariantsWritten(const PgenContext *const pGenContext) {
        return plink2::SpgwGetVidx(pGenContext->spgwp);
    }

    plink2::PgenWriteMode validatePgenWriteMode(const uint32_t pgenWriteModeInt) {
        switch (pgenWriteModeInt) {
            case 0:
            case 1:
            case 2:
                return static_cast<plink2::PgenWriteMode>(pgenWriteModeInt);

            default:
                snprintf(reservedForExceptionMessage,
                         kReservedMessageBufSize,
                         "Invalid pgenWriteMode value (%d), must be one of 0, 1, 2", pgenWriteModeInt);
                throw PgenException(reservedForExceptionMessage);
        }
    }

//cdef class PgenWriter:
//    cdef STPgenWriter* _state_ptr
//    cdef uintptr_t* _nonref_flags
//    cdef PgenGlobalFlags _phase_dosage_gflags
//# preallocate buffers we'll use repeatedly
//    cdef uintptr_t* _genovec
//    cdef uintptr_t* _phasepresent
//    cdef uintptr_t* _phaseinfo
//    cdef uintptr_t* _dosage_present
//    cdef uint16_t* _dosage_main
//
//
//    def __cinit__(self, bytes filename, uint32_t sample_ct,
//                  uint32_t variant_ct, object nonref_flags,
//                  object allele_idx_offsets = None,
//                  bint hardcall_phase_present = False,
//                  bint dosage_present = False,
//                  bint dosage_phase_present = False):
//    if dosage_phase_present and not dosage_present:
//    raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
//    if allele_idx_offsets is not None:
//    for uii in range(variant_ct + 1):
//    if allele_idx_offsets[uii] != uii * 2:
//    raise RuntimeError("Multiallelic variants aren't supported by PgenWriter yet.")
//
//    self._state_ptr = <STPgenWriter*>PyMem_Malloc(sizeof(STPgenWriter))
//    if not self._state_ptr:
//    raise MemoryError()
//    self._nonref_flags = NULL
//    cdef uint32_t nonref_flags_storage = 0
//    cdef uint32_t bitvec_cacheline_ct = DivUp(sample_ct, kBitsPerCacheline)
//    if nonref_flags is not None:
//    if type(nonref_flags) == type(True):
//    if nonref_flags:
//    nonref_flags_storage = 2
//    else:
//    nonref_flags_storage = 1
//    else:
//    nonref_flags_storage = 3
//    if cachealigned_malloc(bitvec_cacheline_ct * kCacheline, &(self._nonref_flags)):
//    raise MemoryError()
//    bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)
//    cdef const char* fname = <const char*>filename
//    cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
//    if hardcall_phase_present:
//    phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
//    if dosage_present:
//    phase_dosage_gflags |= kfPgenGlobalDosagePresent
//            self._phase_dosage_gflags = phase_dosage_gflags
//    assert not dosage_phase_present
//            cdef uintptr_t alloc_cacheline_ct
//    cdef uint32_t max_vrec_len
//            cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, variant_ct, sample_ct, 0, 0, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
//    if reterr != kPglRetSuccess:
//    raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//    cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
//    cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
//    cdef unsigned char* spgw_alloc
//    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
//    raise MemoryError()
//    SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
//    self._genovec = <uintptr_t*>(&(spgw_alloc[alloc_cacheline_ct * kCacheline]))
//    self._phasepresent = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct) * kCacheline]))
//    self._phaseinfo = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + bitvec_cacheline_ct) * kCacheline]))
//    self._dosage_present = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * kCacheline]))
//    self._dosage_main = <uint16_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * kCacheline]))
//    return

//    cpdef append_alleles_batch(self, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_batch, bint all_phased = False):
//    cdef uint32_t batch_size = <uint32_t>allele_int32_batch.shape[0]
//    cdef uintptr_t* genovec = self._genovec
//    cdef int32_t* allele_codes
//            cdef uint32_t uii
//    cdef PglErr reterr
//    if not all_phased:
//    for uii in range(batch_size):
//    allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
//    AlleleCodesToGenoarrUnsafe(allele_codes, NULL, SpgwGetSampleCt(self._state_ptr), genovec, NULL, NULL)
//    reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
//    if reterr != kPglRetSuccess:
//    raise RuntimeError("append_alleles_batch() error " + str(reterr))
//    else:
//    if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//    raise RuntimeError("append_alleles_batch called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
//    for uii in range(batch_size):
//    allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
//    AlleleCodesToGenoarrUnsafe(allele_codes, NULL, SpgwGetSampleCt(self._state_ptr), genovec, NULL, self._phaseinfo)
//    reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, self._phaseinfo, self._state_ptr)
//    if reterr != kPglRetSuccess:
//    raise RuntimeError("append_alleles_batch() error " + str(reterr))
//    return

}