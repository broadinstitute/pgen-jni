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

    //TODO: should we use variant_ct_limit instead of requiring the variant count to be known up front ?
    // Since we need to keep this code in sync with the plink2 pgenlib code, try to use variable names that match
    // the python code to make it easier to track/propagate code changes when updating to newer plink2 commits.
    PgenContext *openPgen(
            const char* cFilename,
            const int pgenWriteModeInt,
            const long numberOfVariants,
            const int sampleCount) {

        // validate the requested pgen write mode
        plink2::PgenWriteMode pgenWriteMode = validatePgenWriteMode(pgenWriteModeInt);

        PgenContext* pGenContext = static_cast<PgenContext *const>(malloc(sizeof(PgenContext)));
        if (pGenContext == nullptr) {
            throw PgenException("Native code failure allocating PgenContext");
        }
        pGenContext->spgwp = static_cast<plink2::STPgenWriter *>(malloc(sizeof(plink2::STPgenWriter)));
        if (pGenContext->spgwp == nullptr) {
            free(pGenContext);
            throw PgenException("Native code failure allocating STPgenWriter");
        }
        //TODO: fix this type
        pGenContext->sampleCount = static_cast<uint32_t>(sampleCount);

        //TODO: fix these types, plink2::DivUp takes two uintptr_t, which are long on Mac
        uint32_t bitvec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kBitsPerCacheline);

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

        //**********
        // TODO: wtf does this do ??
        //          bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)

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
        const plink2::PglErr init1Result = plink2::SpgwInitPhase1(cFilename, //filename
                                                                  nullptr,  // allele index offsets ( for multi allele)
                                                                  nullptr,  // non-ref flags
                                                                  static_cast<uint32_t>(numberOfVariants), // number of variants
                                                                  pGenContext->sampleCount, // sample count
                                                                  0, // optional max allele count
                                                                  pgenWriteMode,
                                                                  plink2::kfPgenGlobal0, //todo- is this right ? type: PgenGlobalFlags phase dosage gflags (genotype?)
                                                                  1, //  non-ref flags storage
                                                                  pGenContext->spgwp, // STPgenWriter * spgwp
                                                                  &alloc_cacheline_ct, //  uintptr_t* alloc_cacheline_ct
                                                                  &max_vrec_len);  // max vrec len ptr
        throwOnPglErr(init1Result, "plink2::SpgwInitPhase1 failed");

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
        //TODO: remove me
        int memSize = (alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + dosage_main_cacheline_ct) *
                      plink2::kCacheline;
        //    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
        if (plink2::cachealigned_malloc(
                (alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) *
                plink2::kCacheline, &spgw_alloc)) {
            throw PgenException("Native code failure allocating plink2::cachealigned_malloc");
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
        int i1 = (pGenContext->sampleCount - 1) / plink2::kBitsPerWordD2;
        pGenContext->genovec[(pGenContext->sampleCount - 1) / plink2::kBitsPerWordD2] = 0; // rely on floor division
        //    self._phasepresent[(sample_ct - 1) // kBitsPerWord] = 0
        int i2 = (pGenContext->sampleCount - 1) / plink2::kBitsPerWord;
        pGenContext->phasepresent[(pGenContext->sampleCount - 1) / plink2::kBitsPerWord] = 0;

        return pGenContext;
    }

    void appendAlleles(const PgenContext *const pGenContext, const int32_t* allele_codes ) {
        //TODO: for now just say allPhased == false;
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
            throw new PgenException("plink2::ConvertMultiAlleleCodesUnsafe return == -1");
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
        throwOnPglErr(pglErr, "appendAlleles multi allelic");
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

        throwOnPglErr(SpgwFinish(pGenContext->spgwp), "Error closing pgen file: SpgwFinish");
        plink2::PglErr cleanupErr;
        plink2::BoolErr bErr = CleanupSpgw(pGenContext->spgwp, &cleanupErr);
        if (bErr) {
            throwOnPglErr(cleanupErr, "Error cleaning up pgen file: CleanupSpgw");
        }
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

// PgenWriter State:
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

//    def __cinit__(self,
//                  bytes filename,
//                  uint32_t sample_ct,
//                  object variant_ct = None,
//                  object nonref_flags = True,
//                  uint32_t allele_ct_limit = 2,
//                  bint hardcall_phase_present = False,
//                  bint dosage_present = False,
//                  bint dosage_phase_present = False,
//                  object variant_ct_limit = None):
//    cdef PgenWriteMode write_mode = kPgenWriteBackwardSeek
//    cdef uint32_t cur_variant_ct_limit
//    if variant_ct is not None:
//      # Could also enforce variant_ct >= variant_ct_limit when both are
//      # provided?
//      cur_variant_ct_limit = variant_ct
//    elif variant_ct_limit is not None:
//      write_mode = kPgenWriteAndCopy
//      cur_variant_ct_limit = variant_ct_limit
//    else:
//      raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct or variant_ct_upper_bound must be provided).")
//    if cur_variant_ct_limit == 0 or cur_variant_ct_limit > kPglMaxVariantCt:
//      raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct must be positive, and less than ~2^31).")
//    if dosage_phase_present and not dosage_present:
//      raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
//    if allele_ct_limit != 2:
//      if (allele_ct_limit < 2) or (allele_ct_limit > kPglMaxAlleleCt):
//          raise RuntimeError("Invalid arguments for PgenWriter constructor (allele_ct_limit must be in [2, 255]).")
//      write_mode = kPgenWriteAndCopy
//
//    self._state_ptr = <STPgenWriter*>PyMem_Malloc(sizeof(STPgenWriter))
//    if not self._state_ptr:
//      raise MemoryError()
//    self._nonref_flags = NULL
//    cdef uint32_t nonref_flags_storage = 0
//    cdef uint32_t bitvec_cacheline_ct = DivUp(sample_ct, kBitsPerCacheline)
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
//          bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)
//    cdef const char* fname = <const char*>filename
//    cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
//    if hardcall_phase_present:
//      phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
//    if dosage_present:
//      phase_dosage_gflags |= kfPgenGlobalDosagePresent
//    self._phase_dosage_gflags = phase_dosage_gflags
//    assert not dosage_phase_present
//
//    cdef uintptr_t alloc_cacheline_ct
//    cdef uint32_t max_vrec_len
//    cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, cur_variant_ct_limit, sample_ct, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
//    if reterr != kPglRetSuccess:
//      raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//    cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
//    cdef uint32_t patch_01_vals_cacheline_ct = DivUp(sample_ct * sizeof(AlleleCode), kCacheline)
//    cdef uint32_t patch_10_vals_cacheline_ct = DivUp(sample_ct * 2 * sizeof(AlleleCode), kCacheline)
//    cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
//    cdef unsigned char* spgw_alloc
//    if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
//      raise MemoryError()
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
//      if all_phased:
//          if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
//              raise RuntimeError("append_alleles called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
//          phaseinfo = self._phaseinfo
//      cdef uint32_t patch_01_ct
//      cdef uint32_t patch_10_ct
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
//      if allele_ct == -1:
//          raise RuntimeError("append_alleles called with invalid allele codes")
//      if <uint32_t>(allele_ct) > allele_ct_limit:
//          raise RuntimeError("append_alleles called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")

//      cdef PglErr reterr
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