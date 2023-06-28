/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

#include "org_broadinstitute_pgen_PgenWriter.h"

#include <iostream>
#include "PgenJniUtils.h"
#include "pgenIO.h"
#include "pgenContext.h"
#include "pgenException.h"
#include "pgenMissingVariantsException.h"
#include "pgenEmptyPgenException.h"

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
                                                 jint writeFlags,
                                                 jlong numberOfVariants,
                                                 jint sampleCount,
                                                 jint maxAltAlleles) {

    // the plink code makes a copy of this filename, so this can be released before this function returns
    const char* const cFilename = env->GetStringUTFChars(filename, nullptr);
 
    jlong pgenHandle;
    try {
        PgenContext* const pgenContext = openPgen(cFilename, pgenWriteModeInt, static_cast<unsigned int>(writeFlags), numberOfVariants, sampleCount, maxAltAlleles);
        env->ReleaseStringUTFChars (filename, cFilename);
        pgenHandle = reinterpret_cast<jlong>(pgenContext);
    } catch (const PgenException& e) {
        env->ReleaseStringUTFChars (filename, cFilename);
        reThrowAsAsyncJavaException(env, e, "Native code failure opening pgen context");
        pgenHandle = 0L;
    }
    return pgenHandle;
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_appendAlleles(JNIEnv *env, jclass object,
                                                      jlong pgenHandle,
                                                      jobject alleleBuffer,
                                                      jobject phaseBuffer,
                                                      jint alleleCount) {
    const int32_t *allele_codes = reinterpret_cast<int32_t*>(env->GetDirectBufferAddress(alleleBuffer));
    if ( !allele_codes ) {
        throwAsyncJavaException(
            env,
            "Native code failure getting address for allele codes in appendAlleles",
            "org/broadinstitute/pgen/PgenException");
        return false;
    } else {
        const unsigned char *phase_buffer = reinterpret_cast<unsigned char*>(env->GetDirectBufferAddress(phaseBuffer));
        if ( !phase_buffer ) {
            throwAsyncJavaException(
                env,
                "Native code failure getting address for phaseBuffer in appendAlleles",
                "org/broadinstitute/pgen/PgenException");
            return false;
        } else {
            PgenContext *pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
            try {
                appendAlleles(pgenContext, allele_codes, phase_buffer, alleleCount);
                return true;
            } catch (const PgenException &e) {
                reThrowAsAsyncJavaException(env, e, "Native code failure in appendAlleles");
                return false;
            }
        }
    }
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_closePgen(JNIEnv *env, jclass object, jlong pgenHandle, jlong droppedVariantCount) {
    PgenContext *pgenContext = reinterpret_cast<PgenContext*>(pgenHandle);
    try {
        closePgen(pgenContext, droppedVariantCount);
        return true;
    } catch (PgenMissingVariantsException &e) {
        // Don't re-throw variant count exceptions as a Java exception, since this function is called from the
        // close method of the Java writer. If the writer was created in a try-with-resources, and writing has
        // terminated prematurely (i.e., in the course of writing the pgen another exception has *already* been
        // thrown), throwing  again from the close method will cause the original exception to be suppressed.
        // So just write the message to stderr and return true.
        std::cerr << "Variant count mismatch detected on close (exception suppressed): " << e.what() << " \n";
        return true;
    } catch (PgenException &e) {
        // Let any other PgenException propagate, but since throwing a Java exception from the close method of
        // the writer can mask a previous exception if it happens in a try-with-resources, log the original
        // error to stderr before we propagate the exception.
        std::cerr << "Error ocurred in  native code during close: " << e.what();
        reThrowAsAsyncJavaException(env, e, "Native code failure closing pgen context");
        throw e;
    } catch (PgenEmptyPgenException &e) {
        // no variants were written - an empty PGEN isn't valid, so give the caller a chance to hande/report that
        throwAsyncJavaException(
            env,
            e.what(),
            "org/broadinstitute/pgen/PgenEmptyPgenException");
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
        throwAsyncJavaException(
            env,
            "Native code failure allocating memory for ByteBuffer",
            "org/broadinstitute/pgen/PgenException");
        return nullptr;
    }
    return env->NewDirectByteBuffer(buf, length);
}

JNIEXPORT jboolean JNICALL
Java_org_broadinstitute_pgen_PgenWriter_destroyByteBuffer( JNIEnv *env, jclass cls, jobject byteBuf ) {
    void *buf = env->GetDirectBufferAddress(byteBuf);
    if ( !buf ) {
        throwAsyncJavaException(
            env,
            "Native code failure getting ByteBuffer address to free",
            "org/broadinstitute/pgen/PgenException");
        return false;
    }
    free(buf);
    return true;
}
