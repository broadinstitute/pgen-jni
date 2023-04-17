#include "org_broadinstitute_pgen_PgenWriter.h"

#include "PgenJniUtils.h"
#include "pgenIO.h"
#include "pgenContext.h"
#include "pgenException.h"
#include "pgenlib_write.h"

using namespace pgenlib;

// Implementation of the JNI access layer. In general this code should do as little as possible,
// only converting to and from Java types, delegating as much as possible to the underlying C++
// pgenlib code.
//
// C++ exceptions from lower layers that are caught here are re-thrown as Java exceptions.

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_pgen_PgenWriter_openPgen (JNIEnv *env, jclass thisObject,
                                                 jstring filename,
                                                 jint pgenWriteModeInt,
                                                 jlong numberOfVariants,
                                                 jint sampleCount) {

    // the plink code makes a copy of filename, so this can be released before this function returns
    const char* const cFilename = env->GetStringUTFChars(filename, nullptr);
 
    jlong pgenHandle;
    try {
        PgenContext* const pgenContext = openPgen(cFilename, pgenWriteModeInt, numberOfVariants, sampleCount);
        env->ReleaseStringUTFChars (filename, cFilename);
        pgenHandle = reinterpret_cast<jlong>(pgenContext);
    } catch (PgenException& e) {
        env->ReleaseStringUTFChars (filename, cFilename);
        reThrowAsJavaException(env, e, "Native code failure opening pgen context");
        pgenHandle = 0L;
    }
    return pgenHandle;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_appendAlleles(JNIEnv* env, jobject object,
                                                      jlong pgenHandle,
                                                      jobject alleleBuffer){
    const int32_t* allele_codes = reinterpret_cast<int32_t*>(env->GetDirectBufferAddress(alleleBuffer));
    if ( !allele_codes ) {
        throwJavaException(env, "Native code failure getting address for seqs ByteBuffer");
    } else {
        PgenContext* pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
        try {
            appendAlleles(pgenContext, allele_codes);
        } catch (PgenException& e) {
            reThrowAsJavaException(env, e, "Native code failure appending alleles");
        }
    }
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_closePgen (JNIEnv * env, jobject object, jlong pgenHandle){
    PgenContext* pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
    try {
        closePgen(pgenContext);
    } catch (PgenException& e) {
        reThrowAsJavaException(env, e, "Native code failure closing pgen context");
    }
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_pgen_PgenWriter_createBuffer( JNIEnv* env, jclass cls, jint length ) {
    void* buf = malloc(length);
    if ( !buf ) {
        throwJavaException(env, "Native code failure allocating memory for  buffer");
        return 0;
    }
    return env->NewDirectByteBuffer(buf, length);
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_destroyByteBuffer( JNIEnv* env, jclass cls, jobject byteBuf ) {
    void* buf = env->GetDirectBufferAddress(byteBuf);
    if ( !buf ) {
        throwJavaException(env, "Native code failure getting ByteBuffer address");
    }
    free(buf);
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
//cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False):
//        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
//        cdef uintptr_t* genovec = self._genovec
//        cdef PglErr reterr
//        if not all_phased:
//            AlleleCodesToGenoarrUnsafe(allele_codes, NULL, SpgwGetSampleCt(self._state_ptr), genovec, NULL, NULL)
//            reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
//        else:
//            AlleleCodesToGenoarrUnsafe(allele_codes, NULL, SpgwGetSampleCt(self._state_ptr), genovec, NULL, self._phaseinfo)
//            reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, self._phaseinfo, self._state_ptr)
//        if reterr != kPglRetSuccess:
//            raise RuntimeError("append_alleles() error " + str(reterr))
//        return
