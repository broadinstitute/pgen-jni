#include <iostream>
#include <cmath>

using namespace std;

#include "pgenContext.h"
#include "pgenException.h"
#include "pgenUtils.h"
#include "pgenlib_misc.h"
#include "pgenlib_write.h"

namespace pgenlib {

    plink2::PgenWriteMode validatePgenWriteMode(const uint32_t anInt);
    static const int kErrMessageBufSize = 1024;

    PgenContext *openPgen(
            const char* cFilename,
            const int pgenWriteModeInt,
            const long variantCount,
            const int sampleCount) {

        // validate the requested pgen write mode, and sample and variant counts
        plink2::PgenWriteMode pgenWriteMode = validatePgenWriteMode(pgenWriteModeInt);
        if (sampleCount < 1) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff, kErrMessageBufSize, "Invalid sample count: %d. At least 1 sample is required.", sampleCount);
            throw PgenException(errMessageBuff);
        } else if (variantCount < 1) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff, kErrMessageBufSize, "Invalid variant count: %ld. At least 1 variant is required.", variantCount);
            throw PgenException(errMessageBuff);
        } else if (variantCount > UINT32_MAX) {
            char errMessageBuff[kErrMessageBufSize];
            snprintf(errMessageBuff, kErrMessageBufSize, "Invalid variant count: %ld exceeds maximum allowable variant count.", variantCount);
            throw PgenException(errMessageBuff);
        }

        PgenContext* pGenContext = static_cast<PgenContext *const>(malloc(sizeof(PgenContext)));
        if (pGenContext == nullptr) {
            throw PgenException("Native code failure allocating PgenContext");
        }
        pGenContext->spgwp = static_cast<plink2::STPgenWriter *>(malloc(sizeof(plink2::STPgenWriter)));
        if (pGenContext->spgwp == nullptr) {
            free(pGenContext);
            throw PgenException("Native code failure allocating STPgenWriter");
        }

        // convert sampleCount and variantCount to the types plink wants
        pGenContext->sampleCount = static_cast<uint32_t>(sampleCount);
        uint32_t variant_ct = static_cast<uint32_t>(variantCount);

        //plink2::DivUp(uintptr_t val, uint32_t divisor)
        uint32_t bitvec_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, plink2::kBitsPerCacheline);

        //    if nonref_flags is not None:
        //      if type(nonref_flags) == type(True):
        //          if nonref_flags:
        //              nonref_flags_storage = 2
        //          else:
        //              nonref_flags_storage = 1
        //      else:
        //          nonref_flags_storage = 3
        //          if cachealigned_malloc(bitvec_cacheline_ct * kCacheline, &(self._nonref_flags)):
        //              raise MemoryError()

        // this code is not using nonref_flags
        // bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)

        //    cdef const char* fname = <const char*>filename
        //    cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
        //    if hardcall_phase_present:
        //      phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
        //    if dosage_present:
        //      phase_dosage_gflags |= kfPgenGlobalDosagePresent
        //    self._phase_dosage_gflags = phase_dosage_gflags
        //    assert not dosage_phase_present

        //    cdef uintptr_t alloc_cacheline_ct
        //    cdef uint32_t max_vrec_len
        //    cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, cur_variant_ct_limit, sample_ct, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
        //    if reterr != kPglRetSuccess:
        //      raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
        uintptr_t alloc_cacheline_ct = 0;
        uint32_t max_vrec_len;
//        PglErr SpgwInitPhase1(
//                const char* __restrict fname,
//                const uintptr_t* __restrict allele_idx_offsets,
//                uintptr_t* __restrict explicit_nonref_flags,
//                uint32_t variant_ct_limit,
//                uint32_t sample_ct,
//                uint32_t allele_ct_upper_bound,
//                PgenWriteMode write_mode,
//                PgenGlobalFlags phase_dosage_gflags,
//                uint32_t nonref_flags_storage,
//                STPgenWriter* spgwp,
//                uintptr_t* alloc_cacheline_ct_ptr,
//                uint32_t* max_vrec_len_ptr)

            const plink2::PglErr init1Result = plink2::SpgwInitPhase1(cFilename, //filename
                                                                  nullptr,  // allele index offsets ( for multi allele)
                                                                  nullptr,  // non-ref flags
                                                                  variant_ct, // number of variants
                                                                  pGenContext->sampleCount, // sample count
                                                                  0, // optional max allele count
                                                                  pgenWriteMode,
                                                                  plink2::kfPgenGlobal0,
                                                                  1, //  non-ref flags storage
                                                                  pGenContext->spgwp, // STPgenWriter * spgwp
                                                                  &alloc_cacheline_ct, //  uintptr_t* alloc_cacheline_ct
                                                                  &max_vrec_len);  // max vrec len ptr
        throwOnPglErr(init1Result, "plink2 initialization (SpgwInitPhase1 failed)");

        //    cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
        //    cdef uint32_t patch_01_vals_cacheline_ct = DivUp(sample_ct * sizeof(AlleleCode), kCacheline)
        //    cdef uint32_t patch_10_vals_cacheline_ct = DivUp(sample_ct * 2 * sizeof(AlleleCode), kCacheline)
        //    cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
        uint32_t genovec_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, plink2::kNypsPerCacheline);
        uint32_t patch_01_vals_cacheline_ct = plink2::DivUp(pGenContext->sampleCount * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t patch_10_vals_cacheline_ct = plink2::DivUp(pGenContext->sampleCount * 2 * sizeof(plink2::AlleleCode), plink2::kCacheline);
        uint32_t dosage_main_cacheline_ct = plink2::DivUp(pGenContext->sampleCount, (2 * plink2::kInt32PerCacheline));

        //    cdef unsigned char* spgw_alloc
        //    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
        //      raise MemoryError()
        //    SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
        // There are two copies of pgenlib.pyx in the plink2 build (that have many differences. One uses + 3, one uses + 5.
        // Prefer the one in src, and go with + 5.
        unsigned char* spgw_alloc;
        //    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
        if (plink2::cachealigned_malloc(
                (alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) *
                plink2::kCacheline, &spgw_alloc)) {
            throw PgenException("Native code failure allocating cachealigned_malloc");
        }
        SpgwInitPhase2(max_vrec_len, pGenContext->spgwp, spgw_alloc);

        //    cdef unsigned char* spgw_alloc_iter = &(spgw_alloc[alloc_cacheline_ct * kCacheline])
        unsigned char* spgw_alloc_iter = &(spgw_alloc[alloc_cacheline_ct * plink2::kCacheline]);
        //    self._allele_ct_limit = allele_ct_limit
        //    self._genovec = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->genovec = (uintptr_t *) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * plink2::kCacheline]);

        //    # Can't skimp on patch_{01,10}_{set,vals} allocations even when
        //    # allele_ct_limit == 2, due to how ConvertMultiAlleleCodesUnsafe()
        //    # works.
        //    # Could skimp on dosage/phase, but that doesn't gain us much.

        //    self._patch_01_set = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->patch_01_set = (uintptr_t*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        //    self._patch_01_vals = <AlleleCode*>(spgw_alloc_iter)
        pGenContext->patch_01_vals = (plink2::AlleleCode*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * plink2::kCacheline]);
        //    self._patch_10_set = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->patch_10_set = (uintptr_t*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        //    self._patch_10_vals = <AlleleCode*>(spgw_alloc_iter)
        pGenContext->patch_10_vals = (plink2::AlleleCode*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * plink2::kCacheline]);
        //    self._phasepresent = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->phasepresent = (uintptr_t *) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        //    self._phaseinfo = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->phaseinfo = (uintptr_t*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        //    self._dosage_present = <uintptr_t*>(spgw_alloc_iter)
        pGenContext->dosage_present = (uintptr_t*) spgw_alloc_iter;
        //    spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
        //    self._dosage_main = <uint16_t*>(spgw_alloc_iter)
        pGenContext->dosage_main = (uint16_t*) spgw_alloc_iter;

        //    # bugfix (16 Apr 2023): SpgwAppendBiallelicGenovec[Hphase] assumes
        //    # trailing bits are clear
        //    self._genovec[(sample_ct - 1) // kBitsPerWordD2] = 0
        pGenContext->genovec[(pGenContext->sampleCount - 1) / plink2::kBitsPerWordD2] = 0; // rely on floor division
        //    self._phasepresent[(sample_ct - 1) // kBitsPerWord] = 0
        pGenContext->phasepresent[(pGenContext->sampleCount - 1) / plink2::kBitsPerWord] = 0;

        return pGenContext;
    }

    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes) {
        //TODO: for now just say allPhased == false since there is no phaseinfo provided;
        bool allPhased = false;
//        uintptr_t* genovec = pGenContext->genovec;
//        plink2::AlleleCodesToGenoarrUnsafe(allele_codes, nullptr, plink2::SpgwGetSampleCt(pGenContext->spgwp), genovec, nullptr, nullptr);
//        plink2::PglErr pglErr = plink2::SpgwAppendBiallelicGenovec(genovec, pGenContext->spgwp);
//        throwOnPglErr(pglErr, "Native code failure adding genotypes");
//    cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False):
//      cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
//      cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
//      cdef uint32_t allele_ct_limit = self._allele_ct_limit
//      cdef uintptr_t* genovec = self._genovec
//      cdef uintptr_t* patch_01_set = self._patch_01_set
//      cdef AlleleCode* patch_01_vals = self._patch_01_vals
//      cdef uintptr_t* patch_10_set = self._patch_10_set
//      cdef AlleleCode* patch_10_vals = self._patch_10_vals
//      cdef uintptr_t* phaseinfo = NULL
        uintptr_t* phaseinfo = NULL;
//      if all_phased:
//          if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//              raise RuntimeError("append_alleles called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
//          phaseinfo = self._phaseinfo
//      cdef uint32_t patch_01_ct
//      cdef uint32_t patch_10_ct
        uint32_t patch_01_ct;
        uint32_t patch_10_ct;
//      cdef int32_t allele_ct = ConvertMultiAlleleCodesUnsafe(
//                                  allele_codes,
//                                  NULL,
//                                  sample_ct,
//                                  genovec,
//                                  patch_01_set,
//                                  patch_01_vals,
//                                  patch_10_set,
//                                  patch_10_vals,
//                                  &patch_01_ct,
//                                  &patch_10_ct,
//                                  NULL,
//                                  phaseinfo)
      int32_t allele_ct = plink2::ConvertMultiAlleleCodesUnsafe(
                                  allele_codes,
                                  NULL,
                                  pGenContext->sampleCount,
                                  pGenContext->genovec,
                                  pGenContext->patch_01_set,
                                  pGenContext->patch_01_vals,
                                  pGenContext->patch_10_set,
                                  pGenContext->patch_10_vals,
                                  &patch_01_ct,
                                  &patch_10_ct,
                                  NULL,
                                  phaseinfo);
//      if allele_ct == -1:
//          raise RuntimeError("append_alleles called with invalid allele codes")
        if (allele_ct == -1) {
            //TODO: would be nice if we could determine what the invalid code is...
            throw PgenException("Attempt to append invalid allele code (plink2::ConvertMultiAlleleCodesUnsafe)");
        }
        //TODO: code this in
//      if <uint32_t>(allele_ct) > allele_ct_limit:
//          raise RuntimeError("append_alleles called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")

//      cdef PglErr reterr
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
                        allele_ct,
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
                      allele_ct,
                      patch_01_ct,
                      patch_10_ct,
                      pGenContext->spgwp);
          }
        }
        throwOnPglErr(pglErr, "appendAlleles");
//      if not all_phased:
//          if (patch_01_ct == 0) and (patch_10_ct == 0):
//              reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
//          else:
//              reterr = SpgwAppendMultiallelicSparse(
//                          genovec,
//                          patch_01_set,
//                          patch_01_vals,
//                          patch_10_set,
//                          patch_10_vals,
//                          allele_ct,
//                          patch_01_ct,
//                          patch_10_ct,
//                          self._state_ptr)
//      else:
//          if (patch_01_ct == 0) and (patch_10_ct == 0):
//              reterr = SpgwAppendBiallelicGenovecHphase(
//                          genovec,
//                          NULL,
//                          phaseinfo,
//                          self._state_ptr)
//          else:
//              reterr = SpgwAppendMultiallelicGenovecHphase(
//                          genovec,
//                          patch_01_set,
//                          patch_01_vals,
//                          patch_10_set,
//                          patch_10_vals,
//                          NULL,
//                          phaseinfo,
//                          allele_ct,
//                          patch_01_ct,
//                          patch_10_ct,
//                          self._state_ptr)
//      if reterr != kPglRetSuccess:
//          raise RuntimeError("append_alleles() error " + str(reterr))
//      return
    }

    void closePgen(const PgenContext *const pGenContext, const long numVariantsDropped) {
        const uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(pGenContext->spgwp);
        const uint32_t writtenVariantCt = plink2::SpgwGetVidx(pGenContext->spgwp);
        const uint32_t droppedVariantCt = static_cast<uint32_t>(numVariantsDropped);

        if ((declaredVariantCt - droppedVariantCt) != writtenVariantCt) {
            //TODO: the plink2 python implementation throws on close if you haven't written as many variants as you
            // initially claimed you would, so we do too (after accounting for variants dropped due to exceeding the
            // maximum allele threshold). but we've come this far - do we REALLY want to throw now ???
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
//    def __cinit__(self, bytes filename, uint32_t sample_ct,
//                  object variant_ct = None, object nonref_flags = True,
//                  uint32_t allele_ct_limit = 2,
//                  bint hardcall_phase_present = False,
//                  bint dosage_present = False,
//                  bint dosage_phase_present = False,
//                  object variant_ct_limit = None):
//        cdef PgenWriteMode write_mode = kPgenWriteBackwardSeek
//        cdef uint32_t cur_variant_ct_limit
//        if variant_ct is not None:
//            # Could also enforce variant_ct >= variant_ct_limit when both are
//            # provided?
//            cur_variant_ct_limit = variant_ct
//        elif variant_ct_limit is not None:
//            write_mode = kPgenWriteAndCopy
//            cur_variant_ct_limit = variant_ct_limit
//        else:
//            raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct or variant_ct_upper_bound must be provided).")
//        if cur_variant_ct_limit == 0 or cur_variant_ct_limit > kPglMaxVariantCt:
//            raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct must be positive, and less than ~2^31).")
//        if dosage_phase_present and not dosage_present:
//            raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
//        if allele_ct_limit != 2:
//            if (allele_ct_limit < 2) or (allele_ct_limit > kPglMaxAlleleCt):
//                raise RuntimeError("Invalid arguments for PgenWriter constructor (allele_ct_limit must be in [2, 255]).")
//            write_mode = kPgenWriteAndCopy
//
//        self._state_ptr = <STPgenWriter*>PyMem_Malloc(sizeof(STPgenWriter))
//        if not self._state_ptr:
//            raise MemoryError()
//        self._nonref_flags = NULL
//        cdef uint32_t nonref_flags_storage = 0
//        cdef uint32_t bitvec_cacheline_ct = DivUp(sample_ct, kBitsPerCacheline)
//        if nonref_flags is not None:
//            if type(nonref_flags) == type(True):
//                if nonref_flags:
//                    nonref_flags_storage = 2
//                else:
//                    nonref_flags_storage = 1
//            else:
//                nonref_flags_storage = 3
//                if cachealigned_malloc(bitvec_cacheline_ct * kCacheline, &(self._nonref_flags)):
//                    raise MemoryError()
//                bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)
//        cdef const char* fname = <const char*>filename
//        cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
//        if hardcall_phase_present:
//            phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
//        if dosage_present:
//            phase_dosage_gflags |= kfPgenGlobalDosagePresent
//        self._phase_dosage_gflags = phase_dosage_gflags
//        assert not dosage_phase_present
//
//        cdef uintptr_t alloc_cacheline_ct
//        cdef uint32_t max_vrec_len
//        cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, cur_variant_ct_limit, sample_ct, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
//        if reterr != kPglRetSuccess:
//            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//        cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
//        cdef uint32_t patch_01_vals_cacheline_ct = DivUp(sample_ct * sizeof(AlleleCode), kCacheline)
//        cdef uint32_t patch_10_vals_cacheline_ct = DivUp(sample_ct * 2 * sizeof(AlleleCode), kCacheline)
//        cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
//        cdef unsigned char* spgw_alloc
//        if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
//            raise MemoryError()
//        SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
//        cdef unsigned char* spgw_alloc_iter = &(spgw_alloc[alloc_cacheline_ct * kCacheline])
//        self._allele_ct_limit = allele_ct_limit
//        self._genovec = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * kCacheline])
//
//        # Can't skimp on patch_{01,10}_{set,vals} allocations even when
//        # allele_ct_limit == 2, due to how ConvertMultiAlleleCodesUnsafe()
//        # works.
//        # Could skimp on dosage/phase, but that doesn't gain us much.
//        self._patch_01_set = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//        self._patch_01_vals = <AlleleCode*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * kCacheline])
//        self._patch_10_set = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//        self._patch_10_vals = <AlleleCode*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * kCacheline])
//        self._phasepresent = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//        self._phaseinfo = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//        self._dosage_present = <uintptr_t*>(spgw_alloc_iter)
//        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
//        self._dosage_main = <uint16_t*>(spgw_alloc_iter)
//        # bugfix (16 Apr 2023): SpgwAppendBiallelicGenovec[Hphase] assumes
//        # trailing bits are clear
//        self._genovec[(sample_ct - 1) // kBitsPerWordD2] = 0
//        self._phasepresent[(sample_ct - 1) // kBitsPerWord] = 0
//        return

//    cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False):
//        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
//        cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
//        cdef uint32_t allele_ct_limit = self._allele_ct_limit
//        cdef uintptr_t* genovec = self._genovec
//        cdef uintptr_t* patch_01_set = self._patch_01_set
//        cdef AlleleCode* patch_01_vals = self._patch_01_vals
//        cdef uintptr_t* patch_10_set = self._patch_10_set
//        cdef AlleleCode* patch_10_vals = self._patch_10_vals
//        cdef uintptr_t* phaseinfo = NULL
//        if all_phased:
//            if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//                raise RuntimeError("append_alleles called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
//            phaseinfo = self._phaseinfo
//        cdef uint32_t patch_01_ct
//        cdef uint32_t patch_10_ct
//        cdef int32_t allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, NULL, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, NULL, phaseinfo)
//        if allele_ct == -1:
//            raise RuntimeError("append_alleles called with invalid allele codes")
//        if <uint32_t>(allele_ct) > allele_ct_limit:
//            raise RuntimeError("append_alleles called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
//        cdef PglErr reterr
//        if not all_phased:
//            if (patch_01_ct == 0) and (patch_10_ct == 0):
//                reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
//            else:
//                reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
//        else:
//            if (patch_01_ct == 0) and (patch_10_ct == 0):
//                reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, phaseinfo, self._state_ptr)
//            else:
//                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, NULL, phaseinfo, allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
//        if reterr != kPglRetSuccess:
//            raise RuntimeError("append_alleles() error " + str(reterr))
//        return

}