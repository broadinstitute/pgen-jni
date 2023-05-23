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
Java_org_broadinstitute_pgen_PgenWriter_openPgen (JNIEnv *env, jclass object,
                                                 jstring filename,
                                                 jint pgenWriteModeInt,
                                                 jlong numberOfVariants,
                                                 jint sampleCount,
                                                 jint maxAltAlleles) {

    // the plink code makes a copy of this filename, so this can be released before this function returns
    const char* const cFilename = env->GetStringUTFChars(filename, nullptr);
 
    jlong pgenHandle;
    try {
        PgenContext* const pgenContext = openPgen(cFilename, pgenWriteModeInt, numberOfVariants, sampleCount, maxAltAlleles);
        env->ReleaseStringUTFChars (filename, cFilename);
        pgenHandle = reinterpret_cast<jlong>(pgenContext);
    } catch (PgenException& e) {
        env->ReleaseStringUTFChars (filename, cFilename);
        reThrowAsJavaException(env, e, "Native code failure opening pgen context");
        pgenHandle = 0L;
    }
    return pgenHandle;
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_appendAlleles(JNIEnv *env, jclass object,
                                                      jlong pgenHandle,
                                                      jobject alleleBuffer,
                                                      jint alleleCount){
    const int32_t *allele_codes = reinterpret_cast<int32_t*>(env->GetDirectBufferAddress(alleleBuffer));
    if ( !allele_codes ) {
        throwJavaException(env, "Native code failure getting address for allele codes in appendAlleles");
    } else {
        PgenContext *pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
        try {
            appendAlleles(pgenContext, allele_codes, alleleCount);
            return true;
        } catch (PgenException &e) {
            reThrowAsJavaException(env, e, "Native code failure in appendAlleles");
            return false;
        }
    }
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_closePgen(JNIEnv *env, jclass object, jlong pgenHandle, jlong droppedVariantCount) {
    PgenContext *pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
    try {
        closePgen(pgenContext, droppedVariantCount);
        return true;
    } catch (PgenException &e) {
        reThrowAsJavaException(env, e, "Native code failure closing pgen context");
        return false;
    }
}

JNIEXPORT jlong JNICALL
Java_org_broadinstitute_pgen_PgenWriter_getPgenVariantCount(JNIEnv *env, jclass object, jlong pgenHandle) {
    PgenContext *pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
    const long varCount = getNumberOfVariantsWritten(pgenContext);
    return varCount;
}

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_pgen_PgenWriter_createBuffer( JNIEnv *env, jclass cls, jint length ) {
    void *buf = malloc(length);
    if ( !buf ) {
        throwJavaException(env, "Native code failure allocating memory for ByteBuffer");
        return nullptr;
    }
    return env->NewDirectByteBuffer(buf, length);
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_destroyByteBuffer( JNIEnv *env, jclass cls, jobject byteBuf ) {
    void *buf = env->GetDirectBufferAddress(byteBuf);
    if ( !buf ) {
        throwJavaException(env, "Native code failure getting ByteBuffer address to free");
        return false;
    }
    free(buf);
    return true;
}
