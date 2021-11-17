#include "org_broadinstitute_pgen_PgenWriter.h"
#include "PgenJniUtils.h"
#include "include/pgenlib_write.h"
#include  "pgenlib_ffi_support.h"
#include <string>

typedef struct BookKeepingStruct {
    plink2::STPgenWriter* spgwp;
    uintptr_t alloc_cacheline_ct_ptr;
    uint32_t max_vrec_len;

    uintptr_t* genovec; //Genotype vector
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    uintptr_t* dosage_present;
    uint16_t* dosage_main;

} BookKeeping;

static void throwErrorMessage( JNIEnv* env, const char* message ) {
    jclass exceptionClass = env->FindClass("org/broadinstitute/pgen/PgenJniException");
    if ( exceptionClass ) env->ThrowNew(exceptionClass, message);
}

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_pgen_PgenWriter_createPgenMetadata (JNIEnv *env, jclass thisObject){
    BookKeeping* bookKeepingPtr = (BookKeeping*)malloc(sizeof(BookKeepingStruct));
    if (bookKeepingPtr == NULL){
       //todo
    }
    bookKeepingPtr->spgwp = (plink2::STPgenWriter*)malloc(sizeof(plink2::STPgenWriter));
    if (bookKeepingPtr->spgwp == NULL){
       //todo
    }
    return (jlong) bookKeepingPtr;
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
Java_org_broadinstitute_pgen_PgenWriter_openPgen (JNIEnv *env, jclass thisObject,
                                                 jstring filename,
                                                 jlong numberOfVariants,
                                                 jlong sampleCount,
                                                 jlong bookKeepingHandle){
    const char* cFilename = env->GetStringUTFChars(filename, NULL);
    BookKeeping* bookKeepingPtr = (BookKeeping*)bookKeepingHandle;

    uintptr_t alloc_cacheline_ct_ptr;
    uint32_t max_vrec_len;

    const plink2::PglErr init1Result = plink2::SpgwInitPhase1(  cFilename, //filename
                     NULL,  // allele index offsets ( for multi allele)
                     NULL,  // non-ref flags
                     (uint32_t) numberOfVariants, // number of variants
                     (uint32_t) sampleCount, // sample count
                      0, // optional max allele count
                      0, // type: PgenGlobalFlags phase dosage gflags (genotype?)
                      0, //  non-ref flags storage
                      bookKeepingPtr->spgwp , // STPgenWriter * spgwp
                      &alloc_cacheline_ct_ptr , //  uintptr_t* alloc_cacheline_ct_ptr
                      &max_vrec_len);  // max vrec len ptr

    uint32_t bitvec_cacheline_ct = plink2::DivUp(sampleCount, plink2::kBitsPerCacheline);

     //        if reterr != kPglRetSuccess:
     //            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
    checkPglErr(env, init1Result, "Initialization phase 1 failed");

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
    SpgwInitPhase2(max_vrec_len, bookKeepingPtr->spgwp, spgw_alloc);
    bookKeepingPtr->genovec = (uintptr_t*)(&(spgw_alloc[alloc_cacheline_ct_ptr * plink2::kCacheline]));
    bookKeepingPtr->phasepresent = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct) * plink2::kCacheline]));
    bookKeepingPtr->phaseinfo = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + bitvec_cacheline_ct) * plink2::kCacheline]));
    bookKeepingPtr->dosage_present = (uintptr_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * plink2::kCacheline]));
    bookKeepingPtr->dosage_main = (uint16_t*)(&(spgw_alloc[(alloc_cacheline_ct_ptr + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * plink2::kCacheline]));
    return 0;
}

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
JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_appendAlleles(JNIEnv* env, jobject object,
                                                      jlong bookKeepingHandle,
                                                      jobject alleleBuffer){
    const int32_t* allele_codes = (int32_t*)env->GetDirectBufferAddress(alleleBuffer);
    if ( !allele_codes ) {
        throwErrorMessage(env, "C code can't get address for seqs ByteBuffer");
        return;
    }
    BookKeeping* bookKeepingPtr = (BookKeeping*)bookKeepingHandle;
    uintptr_t* genovec = bookKeepingPtr->genovec;
    plink2::PglErr reterr;
    plink2::AlleleCodesToGenoarrUnsafe(allele_codes, NULL, plink2::SpgwGetSampleCt(bookKeepingPtr->spgwp), genovec, NULL, NULL);
    reterr = plink2::SpgwAppendBiallelicGenovec(genovec, bookKeepingPtr->spgwp);
    checkPglErr(env, reterr, "Failure while adding genotypes");
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_closePgen (JNIEnv * env, jobject object, jlong bookKeepingHandle){
    BookKeeping* bookKeepingPtr = (BookKeeping*)bookKeepingHandle;
    uint32_t declaredVariantCt = plink2::SpgwGetVariantCt(bookKeepingPtr->spgwp);
    uint32_t writtenVariantCt = plink2::SpgwGetVidx(bookKeepingPtr->spgwp);
    if ( declaredVariantCt != writtenVariantCt) {
        char errbuff[200];
        sprintf(errbuff, "PgenWriter.close() called when number of written variants (%d) unequal to initially declared value (%d)", writtenVariantCt, declaredVariantCt);
        jclass exceptionClass = env->FindClass("org/broadinstitute/pgen/PgenJniException");
        env->ThrowNew(exceptionClass, errbuff);
        return;
    }

   checkPglErr(env, SpgwFinish(bookKeepingPtr->spgwp), "Error while closing PgenFile");
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_pgen_PgenWriter_createBuffer( JNIEnv* env, jclass cls, jint length ) {
    void* buf = malloc(length);
    if ( !buf ) {
        throwErrorMessage(env, "C code can't allocate memory for  buffer");
        return 0;
    }
    return env->NewDirectByteBuffer(buf, length);
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_pgen_PgenWriter_destroyByteBuffer( JNIEnv* env, jclass cls, jobject byteBuf ) {
    void* buf = env->GetDirectBufferAddress(byteBuf);
    if ( !buf ) throwErrorMessage(env, "C code can't get ByteBuffer address");
    free(buf);
}

