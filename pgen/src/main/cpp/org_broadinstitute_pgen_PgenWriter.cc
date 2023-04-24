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
