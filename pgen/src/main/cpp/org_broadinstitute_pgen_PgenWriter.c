#include "org_broadinstitute_pgen_PgenWriter.h"

JNIEXPORT jint JNICALL
Java_org_broadinstitute_pgen_PgenWriter_init(JNIEnv *, jclass, jint value){
    return value;
}

//  def __cinit__(self, bytes filename, uint32_t sample_ct,
//                  uint32_t variant_ct, object nonref_flags,
//                  object allele_idx_offsets = None,
//                  bint hardcall_phase_present = False,
//                  bint dosage_present = False,
//                  bint dosage_phase_present = False):


//        if dosage_phase_present and not dosage_present:
//            raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
//        if allele_idx_offsets is not None:
//            for uii in range(variant_ct + 1):
//                if allele_idx_offsets[uii] != uii * 2:
//                    raise RuntimeError("Multiallelic variants aren't supported by PgenWriter yet.")
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
//        assert not dosage_phase_present
//        cdef uintptr_t alloc_cacheline_ct
//        cdef uint32_t max_vrec_len
//        cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, variant_ct, sample_ct, 0, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)

JNIEXPORT jint JNICALL
Java_org_broadinstitute_pgen_PgenWriter_openPgen (JNIEnv *env, jclass thisObject, jstring filename, jlong numberOfVariants){
    const char* cFilename = env->env->GetStringUTFChars(filename, NULL)
    uint32_t sampleCount = 1
    //             filename  ?      nonrefflags          variantCount         samplecount  ? phase_dosageGflags,
    SpgwInitPhase1(cFilename, NULL, NULL,       (uint32_t) numberOfVariants), sampleCount, 0, )
    return (jint) 0;
}

//        if reterr != kPglRetSuccess:
//            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
//        cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
//        cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
//        cdef unsigned char* spgw_alloc
//        if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
//            raise MemoryError()
//        SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
//        self._genovec = <uintptr_t*>(&(spgw_alloc[alloc_cacheline_ct * kCacheline]))
//        self._phasepresent = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct) * kCacheline]))
//        self._phaseinfo = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + bitvec_cacheline_ct) * kCacheline]))
//        self._dosage_present = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * kCacheline]))
//        self._dosage_main = <uint16_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * kCacheline]))
//        return