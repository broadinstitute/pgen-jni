#ifndef __PGEN_JNI_UTILS__
#define __PGEN_JNI_UTILS__

#include "pgenlib_write.h"
#include "pgenException.h"
#include <jni.h>
#include <iostream>

using namespace pgenlib;

static int ERROR_MESSAGE_BUFFER_SIZE = 1024;

// Throw a Java exception (PgenJniException) with the given error message. Note that control RETURNS
// to the caller after the exception is thrown.
static bool throwJavaException( JNIEnv* env, const char* message ) {
    jclass exceptionClass = env->FindClass("org/broadinstitute/pgen/PgenJniException");
    bool result = false;

    if ( exceptionClass ) {
        jint throwResult = env->ThrowNew(exceptionClass, message);
        if (throwResult < 0) {
            std::cerr << "Failure throwing Java exception from native code while handling underlying exception caused by: " << message;
        } else {
            result = true;
        }
    } else {
        std::cerr << "Unable to find Java exception class while handling underlying exception caused by: " << message;
    }
    return result;
}

// Re-throw a PgenException (that originated in underlying *PGENLIB* code) as a Java exception.
//  Note that control RETURNS to the caller after the exception is thrown.
static bool reThrowAsJavaException( JNIEnv* env, PgenException& pgenException, const char* context) {
    //TODO: construct a proper error message using context:
    //return throwErrorMessage(env, context + " : " + pgenException.what());
    return throwJavaException(env, pgenException.what());
}

#endif // __PGEN_JNI_UTILS__