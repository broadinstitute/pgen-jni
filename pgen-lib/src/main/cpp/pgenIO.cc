#include <cmath>

#include "pgenContext.h"
#include "pgenException.h"
#include "pgenMissingVariantsException.h"
#include "pgenEmptyPgenException.h"
#include "pgenUtils.h"
#include "pgenIO.h"
#include "pgenlib_misc.h"
#include "pgenlib_write.h"

namespace pgenlib {
    static const int kErrMessageBufSize = 1024;
    static plink2::PgenWriteMode ValidatePgenWriteMode(const uint32_t pGenWriteMode, const long variantCount);
    static plink2::PgenGlobalFlags PgenlibFlagsToPlink2Flags(const uint32_t pgenlibFlags);
    static PgenContext *InitPgenContext(
        const char* cFilename,
        const plink2::PgenWriteMode pgenWriteMode,
        const uint32_t writeFlags,
        const long variantCount,
        const int sampleCount,
        const int maxAltAlleles);
    static bool GetAllPhased(const PgenContext *pGenContext, const unsigned char *phase_bytes);
    static void AppendAllelesPartiallyPhased(
            const PgenContext *const pGenContext,
            const int32_t* allele_codes,
            const unsigned char* phase_bytes,
            const int32_t allele_ct);
    static void AppendAllelesAllOrNonePhased(
            const PgenContext *const pGenContext,
            const int32_t* allele_codes,
            const unsigned char* phase_bytes,
            const int32_t allele_ct,
            const bool allPhased);

    /**
     * Start a new PGEN write session, and return a pointer to a PgenContext for the writer.
     *
     * The PgenContext can be used to write an entire pgen file. depending on the PGEN file mode used (see below)
     * a .pgen.pgi file may also be created.
     *
     * Only diploid genomes are supported.
     *
     * The allele codes for the genotypes for each variant can be provided by passing the PgenContext returned by
     * this function to a series of calls to appendAlleles, after which the PgenContext session should be closed
     * via a call to closePgen.
     *
     * An example PGEN writer lifecycle is illustrated here:
     *
     *      const pgenlib::PgenContext *const pgen_context = pgenlib::openPgen(
     *          file_name,
     *          pgen_write_mode,
     *          pgen_write_flags,
     *          n_variants,
     *          n_samples,
     *          plink2::kPglMaxAltAlleleCt);
     *
     *      for (int i = 0; i < n_variants; i++) {
     *          pgenlib::AppendAlleles(pgen_context, allele_codes, allele_ct);
     *      }
     *      long variantCount = GetNumberOfVariantsWritten(pgen_context);
     *      ClosePgen(pgen_context, 0);
     *
     *  Once the PgenContext has been closed, it can no longer be used to write allele codes.
     *
     * @param cFilename - the pgen file to write
     * @param pgenWriteModeInt - unsigned integer representing the file mode, with permitted values drawn from integer
     * values of plink2::PgenWriteMode (1, 2 or 3). An exception will be thrown if any other value is provided. this
     * determines the pgen file mode that is used (i.e, whether there is a separate .pgi index)
     * @param writeFlags - unsigned integer bitwise write flags, with valid values drawn from {kWriteFlagPreservePhasing,
     * kWriteFlagMultiAllelic}. kWriteFlagPreservePhasing should only be used if phasing information is present in the
     * source genotypes and a phasing track must be provided when calling appendAlleles. kWriteFlagMultiAllelic should
     * be included if multi-allelic genotypes are present. kWriteFlagMultiAllelic should only be used when
     * kWriteFlagPreservePhasing is used (!).
     * @param variantCount - the number of variants to be written. if fewer variants are written, an exception will
     * be thrown when the writer is closed by a call to closePgen. must be in the range 1..plink2::kPglMaxVariantCt
     * @param sampleCount - the number of samples (genotypes) in the data set. Must be > 0.
     * @param maxAltAlleles - the maximum number of alleles for any variant that will be written - this determines the
     * range of valid (zero based) allele codes that can be provided when writing genotypes to this writer. Must be
     * in the range 2..plink2::kPglMaxAltAlleleCt. If the variant count is unknown when the writer is created, use
     * the value pgenlib::kVariantCountUnknown (although in this case, write mode
     * plink2::PgenWriteMode::kPgenWriteBackwardSeek (3) may not be used).
     *
     * @return a PgenContext
     */
    PgenContext *OpenPgen(
            const char* cFilename,
            const uint32_t pgenWriteModeInt,
            const uint32_t writeFlags,
            const long variantCount,
            const int sampleCount,
            const int maxAltAlleles) {

        // validate the requested pgen write mode, and sample and variant counts
        plink2::PgenWriteMode pgenWriteMode = ValidatePgenWriteMode(pgenWriteModeInt, variantCount);
        if (sampleCount < 1) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid sample count: %d. At least 1 sample is required.",
                     sampleCount);
            throw PgenException(errMessageBuff); // PgenException makes a copy of errMessageBuff
        } else if (variantCount < 1L) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff,
                     kErrMessageBufSize,
                     "Invalid variant count: %ld. Variant count must be > 0.",
                     variantCount);
            throw PgenException(errMessageBuff);  // PgenException makes a copy of errMessageBuff
        } else if (variantCount > static_cast<long>(plink2::kPglMaxVariantCt)) {
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

        // TODO: this seems weird, but according to the comment in ...., only set multiallelic flag if there is
        // also phasing info, even if you really have multi-allelics ? So, if you have multi-allelic data but
        // no phasing, don't set the multi-allelic bit ? really ?
        // Should we relax this weirdness in the API, and instead just always accept the kWriteFlagMultiAllelic, but
        // when it's present, silently remove it before delegating to plink in the C code if kWriteFlagPreservePhasing
        // isn't also set ?
        if ((writeFlags & kWriteFlagMultiAllelic) && !(writeFlags & kWriteFlagPreservePhasing)) {
            throw PgenException("The multi-allelic write flag should only be used if phasing information is also provided (even if the underlying data is multiallelic).");
        }

        return InitPgenContext(cFilename, pgenWriteMode, writeFlags, variantCount, sampleCount, maxAltAlleles);
    }

    PgenContext *InitPgenContext(
            const char* cFilename,
            const plink2::PgenWriteMode pgenWriteMode,
            const uint32_t writeFlags,
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
        pGenContext->write_flags = writeFlags;

        // convert sampleCount and variantCount to the types plink uses
        pGenContext->sample_count = static_cast<uint32_t>(sampleCount);
        uint32_t variant_ct = static_cast<uint32_t>(variantCount);
        // total max allele count is max alt count + 1
        pGenContext->allele_ct_limit = static_cast<uint32_t>(maxAltAlleles + 1);

        uint32_t bitvec_cacheline_ct = plink2::DivUp(pGenContext->sample_count, plink2::kBitsPerCacheline);
        uintptr_t alloc_cacheline_ct = 0;
        const plink2::PglErr init1Result = plink2::SpgwInitPhase1(cFilename,
                                                                nullptr,  // allele index offsets (for reading multi allele ?)
                                                                nullptr, // non-ref flags
                                                                variant_ct,
                                                                pGenContext->sample_count,
                                                                pGenContext->allele_ct_limit,
                                                                pgenWriteMode,
                                                                PgenlibFlagsToPlink2Flags(pGenContext->write_flags),
                                                                1, // non-ref flags storage
                                                                pGenContext->spgwp,
                                                                &alloc_cacheline_ct,
                                                                &pGenContext->max_vrec_len);
        throwOnPglErr(init1Result, "plink2 initialization (SpgwInitPhase1 failed)");

        uint32_t genovec_cacheline_ct = plink2::DivUp(pGenContext->sample_count, plink2::kNypsPerCacheline);
        uint32_t patch_01_vals_cacheline_ct = plink2::DivUp(pGenContext->sample_count * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t patch_10_vals_cacheline_ct = plink2::DivUp(pGenContext->sample_count * 2 * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t dosage_main_cacheline_ct = plink2::DivUp(pGenContext->sample_count, (2 * plink2::kInt32PerCacheline));

        // There are two copies of pgenlib.pyx in the plink2 build, and they have many differences. One uses +3 for
        // this calculation, and one uses +5. Prefer the one in src (since thats the one that is the template for this
        // code), and go with +5.
        // Keep the pointer to the arena block in pGenContext so we can free it at the end.
        if (plink2::cachealigned_malloc(
                (alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct +
                patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) *
                plink2::kCacheline, &pGenContext->spgw_alloc)) {
            throw PgenException("Native code failure (cachealigned_malloc) allocating spgw_alloc");
        }
        SpgwInitPhase2(pGenContext->max_vrec_len, pGenContext->spgwp, pGenContext->spgw_alloc);

        unsigned char* spgw_alloc_iter = &(pGenContext->spgw_alloc[alloc_cacheline_ct * plink2::kCacheline]);

        pGenContext->genovec = (uintptr_t *) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * plink2::kCacheline]);
        // # Can't skimp on patch_{01,10}_{set,vals} allocations even when
        // # allele_ct_limit == 2, due to how ConvertMultiAlleleCodesUnsafe()
        // # works.
        // # Could skimp on dosage/phase, but that doesn't gain us much.
        pGenContext->patch_01_set = (uintptr_t*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_01_vals = (plink2::AlleleCode*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_10_set = (uintptr_t*) spgw_alloc_iter;
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        pGenContext->patch_10_vals = (plink2::AlleleCode*) spgw_alloc_iter;

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
        pGenContext->genovec[(pGenContext->sample_count - 1) / plink2::kBitsPerWordD2] = 0; // floor division
        pGenContext->phasepresent[(pGenContext->sample_count - 1) / plink2::kBitsPerWord] = 0; // floor division

        return pGenContext;
    }

    // Since this code is ported from the plink2 python code, we try to retain the same structure as that
    // code in order to make it  easier to find and integrate changes as they are made to plink2. However, currently
    // the python code doesn't seem to be quite right, specifically, although there appears to be a single unified
    // code path for handing non-phased, mixed phase and all-phased data that works (in the generic "appendAlleles"
    // function), that function doesn't seem to work in the case where you have mixed phase. Until that is resolved,
    // we have to retain two very similar, but slightly different functions from python ("appendAlleles",
    // named "appendAllelesNotPartiallyPhased" here, and "appendAllelesPartiallyPhased", named
    // "appendAllelesPartiallyPhased") here. The top-level "appendAlleles" function decides which of those to call.
    /**
     * Append one variant's worth of allele code (genotypes) to a pgen file.
     * @param pGenContext - the PgenContext for the writer
     * @param allele_codes - array of allele codes to be written
     * @param phase_bytes - phasing (1 for phased, 0 for not phased). must be present when kWriteFlagPreservePhasing was
     * used to create the PgenWriter, otherwise ignored (may be null)
     * @param allele_ct - the number of possible allele values for this variant (not the number of unique alleles
     * that are ACTUALLY observed/present in allele_codes)
   */
   void AppendAlleles(
           const PgenContext *const pGenContext,
           const int32_t* allele_codes,
           const unsigned char* phase_bytes,
           const int32_t allele_ct) {

         // determine up front whether all the genotypes are phased so we can take the right code path through plink
         //TODO: we could probably skip this pass through the phasing track altogether if we required the caller to
         // keep track of the "allPhased" state while assembling the phasing data, and then provide it via a parameter
         bool allPhased = GetAllPhased(pGenContext, phase_bytes);
         if (pGenContext->write_flags & kWriteFlagPreservePhasing && !allPhased) {
             // there is a phasing track, but the genotypes are either mixed phase or all un-phased
             AppendAllelesPartiallyPhased(
                     pGenContext,
                     allele_codes,
                     phase_bytes,
                     allele_ct);
         } else {
             // there is either no phasing track, or there is a phasing track and all the genotypes are phased
             AppendAllelesAllOrNonePhased(
                     pGenContext,
                     allele_codes,
                     phase_bytes,
                     allele_ct,
                     allPhased);
         }
     }

    void AppendAllelesPartiallyPhased(
            const PgenContext *const pGenContext,
            const int32_t* allele_codes,
            const unsigned char* phase_bytes,
            const int32_t allele_ct) {
        uint32_t patch_01_ct;
        uint32_t patch_10_ct;
        int32_t observed_allele_ct = plink2::ConvertMultiAlleleCodesUnsafe(
                allele_codes,
                phase_bytes, //  may be null
                pGenContext->sample_count,
                pGenContext->genovec,
                pGenContext->patch_01_set,
                pGenContext->patch_01_vals,
                pGenContext->patch_10_set,
                pGenContext->patch_10_vals,
                &patch_01_ct,
                &patch_10_ct,
                pGenContext->phasepresent,
                pGenContext->phaseinfo);
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

        plink2::PglErr pglErr;
        if ((patch_01_ct == 0) and (patch_10_ct == 0)) {
            pglErr = SpgwAppendBiallelicGenovecHphase(
                    pGenContext->genovec,
                    pGenContext->phasepresent,
                    pGenContext->phaseinfo,
                    pGenContext->spgwp);
        } else {
            pglErr = SpgwAppendMultiallelicGenovecHphase(
                    pGenContext->genovec,
                    pGenContext->patch_01_set,
                    pGenContext->patch_01_vals,
                    pGenContext->patch_10_set,
                    pGenContext->patch_10_vals,
                    pGenContext->phasepresent,
                    pGenContext->phaseinfo,
                    write_allele_ct,
                    patch_01_ct,
                    patch_10_ct,
                    pGenContext->spgwp);
        }
        throwOnPglErr(pglErr, "appendAlleles");
    }

   void AppendAllelesAllOrNonePhased(
             const PgenContext *const pGenContext,
             const int32_t* allele_codes,
             const unsigned char* phase_bytes,
             const int32_t allele_ct,
             const bool allPhased) {
        uint32_t patch_01_ct;
        uint32_t patch_10_ct;
        int32_t observed_allele_ct = plink2::ConvertMultiAlleleCodesUnsafe(
                allele_codes,
                phase_bytes, //  may be null
                pGenContext->sample_count,
                pGenContext->genovec,
                pGenContext->patch_01_set,
                pGenContext->patch_01_vals,
                pGenContext->patch_10_set,
                pGenContext->patch_10_vals,
                &patch_01_ct,
                &patch_10_ct,
                pGenContext->phasepresent,
                pGenContext->phaseinfo);
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

        plink2::PglErr pglErr;
        if (!allPhased) {
            if ((patch_01_ct == 0) and (patch_10_ct == 0)) {
                pglErr = SpgwAppendBiallelicGenovec(pGenContext->genovec, pGenContext->spgwp);
            } else {
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
                      pGenContext->phasepresent,
                      pGenContext->phaseinfo,
                      pGenContext->spgwp);
          } else {
              pglErr = SpgwAppendMultiallelicGenovecHphase(
                      pGenContext->genovec,
                      pGenContext->patch_01_set,
                      pGenContext->patch_01_vals,
                      pGenContext->patch_10_set,
                      pGenContext->patch_10_vals,
                      pGenContext->phasepresent,
                      pGenContext->phaseinfo,
                      write_allele_ct,
                      patch_01_ct,
                      patch_10_ct,
                      pGenContext->spgwp);
          }
        }
        throwOnPglErr(pglErr, "appendAlleles");
    }

    /**
     * Close a pGenContext, flush the output, and close the pgen file. The pGenContext is no longer valid after this call.
     *
     * @param pGenContext - the pgen context for this writer
     * @param numVariantsDropped - the number of variants dropped (the number f variants dropped, plus the number
     * of variants written, must equal the number of variants projected to be written when the pgen context was
     * initially opened. otherwise a pgenlib::PGenException will be thrown.
     */
    void ClosePgen(const PgenContext *const pGenContext, const long numVariantsDropped) {
        const uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(pGenContext->spgwp);
        const uint32_t writtenVariantCt = plink2::SpgwGetVidx(pGenContext->spgwp);
        const uint32_t droppedVariantCt = static_cast<uint32_t>(numVariantsDropped);

        if ((declaredVariantCt != static_cast<long>(pgenlib::kVariantCountUnknown)) &&
            ((declaredVariantCt - droppedVariantCt) != writtenVariantCt)) {
            // the plink2 python implementation of the pgen writer throws on close if you haven't written as many
            // variants as you initially claimed you would (at least in the case where the variant count is known
            // up front), so we do the same here (after accounting for variants dropped due to exceeding the
            // maximum allele threshold).
            char errMessage[kErrMessageBufSize];
            snprintf(errMessage,
                     kErrMessageBufSize,
                      "closePgen called with number of variants written (%d) not equal to (declared - dropped)  (%d - %d) = %d",
                      writtenVariantCt,
                      declaredVariantCt,
                      droppedVariantCt,
                      declaredVariantCt - droppedVariantCt);
            // throw a C++ exception that is specific to this (variant count) failure, so the Java code for the
            // writer can catch and handle that case without propagating it, because throwing from the Closeable
            // "close" method, when the writer is created within a try-with-resources statement can mask other
            // exceptions
            throw PgenMissingVariantsException(errMessage);
        } else if (writtenVariantCt != 0) {
            // guard against calling the plink2 finish/cleanup methods in the case where no writes have been made
            // because doing so triggers asserts in the plink code, presumably because downstream code paths can't
            // handle it:
            // Assertion failed: (variant_ct), function PwcFinish, file pgenlib_write.cc, line 2284.
            throwOnPglErr(SpgwFinish(pGenContext->spgwp), "Error closing pgen file: SpgwFinish");

            // there may be a bug in plink2 pgen-lib, since I think the plink2 VCF importer only does one or the
            // other of SpgwFinish and CleanupSpgw (SpgwFinish on success, CleanupSpgw on failure), but not both.
            // But if we don't do both here, the output doesn't seem to get flushed until the process exits.
            plink2::PglErr cleanupErr;
            plink2::BoolErr bErr = CleanupSpgw(pGenContext->spgwp, &cleanupErr);
            if (bErr) {
                throwOnPglErr(cleanupErr, "Error cleaning up on pgen close: CleanupSpgw");
            }
        }

        free(pGenContext->spgwp);
        plink2::aligned_free(pGenContext->spgw_alloc);
        free(reinterpret_cast<void *>(const_cast<PgenContext *>(pGenContext)));

        if (writtenVariantCt == 0) {
            throw PgenEmptyPgenException(
                    "An empty PGEN is not valid - at least one variant site must be written to a PGEN. The PGEN file is not valid");
        }
    }

    long GetNumberOfVariantsWritten(const PgenContext *const pGenContext) {
        return plink2::SpgwGetVidx(pGenContext->spgwp);
    }

    plink2::PgenWriteMode ValidatePgenWriteMode(const uint32_t pgenWriteModeInt, const long variantCount) {
        switch (pgenWriteModeInt) {
            case plink2::kPgenWriteBackwardSeek:
                if (variantCount == static_cast<long>(pgenlib::kVariantCountUnknown)) {
                    char errMessageBuff[kErrMessageBufSize];
                    snprintf(errMessageBuff,
                             kReservedMessageBufSize,
                             "pgenWriteMode value (%u) requires a known variant count, and cannot be used with the unknown variant count sentinel value (%d)",
                             pgenWriteModeInt,
                             plink2::kPglMaxVariantCt);
                    throw PgenException(errMessageBuff);
                }
                // fall through
            case plink2::kPgenWriteSeparateIndex:
            case plink2::kPgenWriteAndCopy:
                return static_cast<plink2::PgenWriteMode>(pgenWriteModeInt);

            default:
                snprintf(reservedForExceptionMessage,
                         kReservedMessageBufSize,
                         "Invalid pgenWriteMode value (%u), must be one of 0, 1, 2", pgenWriteModeInt);
                throw PgenException(reservedForExceptionMessage);
        }
    }

    // assumes that the pgenlibFlags have been properly validated, and that kWriteFlagMultiAllelic is only set if
    // there is also phasing info
    plink2::PgenGlobalFlags PgenlibFlagsToPlink2Flags(const uint32_t pgenlibFlags) {
        plink2::PgenGlobalFlags plinkFlags = plink2::kfPgenGlobal0;
        if (pgenlibFlags & pgenlib::kWriteFlagMultiAllelic) {
            plinkFlags |= plink2::kfPgenGlobalMultiallelicHardcallFound;
        }
        if (pgenlibFlags & pgenlib::kWriteFlagPreservePhasing) {
            plinkFlags |= plink2::kfPgenGlobalHardcallPhasePresent;
        }
        return plinkFlags;
    }

    // In order to take the correct code path through plink, we need to determine up front whether all
    // the genotypes are phased. If there is no phasing track at all, allPhased = false, but otherwise,
    // start out assuming all are phased and then bail as soon as we find one unphased genotype.
    bool GetAllPhased(const PgenContext *pGenContext, const unsigned char *phase_bytes) {
        bool allPhased = false;
        if (pGenContext->write_flags & kWriteFlagPreservePhasing) {
            if (phase_bytes == nullptr) {
                throw PgenException("A phasing track is required since kWriteFlagPreservePhasing was specified");
            }
            allPhased = true;
            for (int i = 0; i < pGenContext->sample_count;  i++) {
                if (!phase_bytes[i]) {
                    return false;
                }
            }
        }
        return allPhased;
    }

/***********************************************************************************************************
 * The Python source below is the template for the C++ implementation in this file, and is taken from plink2
 * file "2.0/Python/src/pgenlib/pgenlib.pyx" in the plink2 repo https://github.com/chrchang/plink-ng/. NOTE
 * that there is more than one copy/version of the file named "pgenlib.pyx" in that repo, and the one in
 * "2.0/Python/src/pgenlib/pgenlib.pyx" is the once used here, since it appears to be the most up to date
 * at this time.
 */

//cdef class PgenWriter:
//    cdef STPgenWriter* _state_ptr
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
//
//
//def __cinit__(self, bytes filename, uint32_t sample_ct,
//              object variant_ct = None, object nonref_flags = True,
//              uint32_t allele_ct_limit = 2,
//              bint hardcall_phase_present = False,
//              bint dosage_present = False,
//              bint dosage_phase_present = False,
//              object variant_ct_limit = None):
//    cdef PgenWriteMode write_mode = kPgenWriteBackwardSeek
//    cdef uint32_t cur_variant_ct_limit
//    if variant_ct is not None:
//        # Could also enforce variant_ct >= variant_ct_limit when both are
//        # provided?
//        cur_variant_ct_limit = variant_ct
//    elif variant_ct_limit is not None:
//        write_mode = kPgenWriteAndCopy
//        cur_variant_ct_limit = variant_ct_limit
//    else:
//        raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct or variant_ct_upper_bound must be provided).")
//    if cur_variant_ct_limit == 0 or cur_variant_ct_limit > kPglMaxVariantCt:
//        raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct must be positive, and less than ~2^31).")
//    if dosage_phase_present and not dosage_present:
//        raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
//    if allele_ct_limit != 2:
//        if (allele_ct_limit < 2) or (allele_ct_limit > kPglMaxAlleleCt):
//            raise RuntimeError("Invalid arguments for PgenWriter constructor (allele_ct_limit must be in [2, 255]).")
//        write_mode = kPgenWriteAndCopy
//
//    self._state_ptr = <STPgenWriter*>PyMem_Malloc(sizeof(STPgenWriter))
//    if not self._state_ptr:
//        raise MemoryError()
//    self._nonref_flags = NULL
//    cdef uint32_t nonref_flags_storage = 0
//    cdef uint32_t bitvec_cacheline_ct = DivUp(sample_ct, kBitsPerCacheline)
//    if nonref_flags is not None:
//        if type(nonref_flags) == type(True):
//            if nonref_flags:
//                nonref_flags_storage = 2
//            else:
//                nonref_flags_storage = 1
//        else:
//            nonref_flags_storage = 3
//            if cachealigned_malloc(bitvec_cacheline_ct * kCacheline, &(self._nonref_flags)):
//                raise MemoryError()
//            bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)
//    cdef const char* fname = <const char*>filename
//    cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
//    if hardcall_phase_present:
//        phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
//    if dosage_present:
//        phase_dosage_gflags |= kfPgenGlobalDosagePresent
//    self._phase_dosage_gflags = phase_dosage_gflags
//    assert not dosage_phase_present
//
//    cdef uintptr_t alloc_cacheline_ct
//    cdef uint32_t max_vrec_len
//    cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, cur_variant_ct_limit, sample_ct, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
//    if reterr != kPglRetSuccess:
//        raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//    cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
//    cdef uint32_t patch_01_vals_cacheline_ct = DivUp(sample_ct * sizeof(AlleleCode), kCacheline)
//    cdef uint32_t patch_10_vals_cacheline_ct = DivUp(sample_ct * 2 * sizeof(AlleleCode), kCacheline)
//    cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
//    cdef unsigned char* spgw_alloc
//    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
//        raise MemoryError()
//    SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
//    cdef unsigned char* spgw_alloc_iter = &(spgw_alloc[alloc_cacheline_ct * kCacheline])
//    self._allele_ct_limit = allele_ct_limit
//    self._genovec = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * kCacheline])
//
//    # Can't skimp on patch_{01,10}_{set,vals} allocations even when
//    # allele_ct_limit == 2, due to how ConvertMultiAlleleCodesUnsafe()
//    # works.
//    # Could skimp on dosage/phase, but that doesn't gain us much.
//    self._patch_01_set = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//    self._patch_01_vals = <AlleleCode*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * kCacheline])
//    self._patch_10_set = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//    self._patch_10_vals = <AlleleCode*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * kCacheline])
//    self._phasepresent = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//    self._phaseinfo = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//    self._dosage_present = <uintptr_t*>(spgw_alloc_iter)
//    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//    self._dosage_main = <uint16_t*>(spgw_alloc_iter)
//    # bugfix (16 Apr 2023): SpgwAppendBiallelicGenovec[Hphase] assumes
//    # trailing bits are clear
//    self._genovec[(sample_ct - 1) // kBitsPerWordD2] = 0
//    self._phasepresent[(sample_ct - 1) // kBitsPerWord] = 0
//    return

//cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False, object allele_ct = None):
//    cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
//    cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
//    cdef uint32_t allele_ct_limit = self._allele_ct_limit
//    cdef uintptr_t* genovec = self._genovec
//    cdef uintptr_t* patch_01_set = self._patch_01_set
//    cdef AlleleCode* patch_01_vals = self._patch_01_vals
//    cdef uintptr_t* patch_10_set = self._patch_10_set
//    cdef AlleleCode* patch_10_vals = self._patch_10_vals
//    cdef uintptr_t* phaseinfo = NULL
//    if all_phased:
//        if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//            raise RuntimeError("append_alleles called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
//        phaseinfo = self._phaseinfo
//    cdef uint32_t patch_01_ct
//    cdef uint32_t patch_10_ct
//    cdef int32_t observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, NULL, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, NULL, phaseinfo)
//    if observed_allele_ct == -1:
//        raise RuntimeError("append_alleles called with invalid allele codes")
//    cdef uint32_t write_allele_ct = <uint32_t>(observed_allele_ct)
//    if write_allele_ct > allele_ct_limit:
//        raise RuntimeError("append_alleles called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
//    if allele_ct is not None:
//        if allele_ct < write_allele_ct:
//            raise RuntimeError("append_alleles called with allele codes >= allele_ct argument")
//        if allele_ct > allele_ct_limit:
//            raise RuntimeError("append_alleles called with allele_ct > allele_ct_limit")
//        write_allele_ct = allele_ct
//    cdef PglErr reterr
//    if not all_phased:
//        if (patch_01_ct == 0) and (patch_10_ct == 0):
//            reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
//        else:
//            reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
//    else:
//        if (patch_01_ct == 0) and (patch_10_ct == 0):
//            reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, phaseinfo, self._state_ptr)
//        else:
//            reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, NULL, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
//    if reterr != kPglRetSuccess:
//        raise RuntimeError("append_alleles() error " + str(reterr))
//    return
//
//cpdef append_partially_phased(self, np.ndarray[np.int32_t,mode="c"] allele_int32, np.ndarray[np.uint8_t,cast=True] phasepresent, object allele_ct = None):
//    if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//        raise RuntimeError("append_partially_phased cannot be called when PgenWriter was constructed with hardcall_phase_present False")
//    cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
//    cdef unsigned char* phasepresent_bytes = <unsigned char*>(&(phasepresent[0]))
//    cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
//    cdef uint32_t allele_ct_limit = self._allele_ct_limit
//    cdef uintptr_t* genovec = self._genovec
//    cdef uintptr_t* phasepresent_buf = self._phasepresent
//    cdef uintptr_t* phaseinfo = self._phaseinfo
//    cdef uintptr_t* patch_01_set = self._patch_01_set
//    cdef AlleleCode* patch_01_vals = self._patch_01_vals
//    cdef uintptr_t* patch_10_set = self._patch_10_set
//    cdef AlleleCode* patch_10_vals = self._patch_10_vals
//    cdef uint32_t patch_01_ct
//    cdef uint32_t patch_10_ct
//    cdef int32_t observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, phasepresent_bytes, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, phasepresent_buf, phaseinfo)
//    if observed_allele_ct == -1:
//        raise RuntimeError("append_partially_phased called with invalid allele codes")
//    cdef uint32_t write_allele_ct = <uint32_t>(observed_allele_ct)
//    if write_allele_ct > allele_ct_limit:
//        raise RuntimeError("append_partially_phased called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
//    if allele_ct is not None:
//        if allele_ct < write_allele_ct:
//            raise RuntimeError("append_partially_phased called with allele codes >= allele_ct argument")
//        if allele_ct > allele_ct_limit:
//            raise RuntimeError("append_partially_phased called with allele_ct > allele_ct_limit")
//        write_allele_ct = allele_ct
//    cdef PglErr reterr
//    if (patch_01_ct == 0) and (patch_10_ct == 0):
//        reterr = SpgwAppendBiallelicGenovecHphase(genovec, phasepresent_buf, phaseinfo, self._state_ptr)
//    else:
//        reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent_buf, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
//    if reterr != kPglRetSuccess:
//        raise RuntimeError("append_partially_phased() error " + str(reterr))
//    return

}