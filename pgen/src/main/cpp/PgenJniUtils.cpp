#include "PgenJniUtils.h"
#include <iostream>

using namespace pgenlib;
using namespace std;

// Throw a Java exception (PgenJniException) with the given error message. Note that control RETURNS
// to the caller after the exception is thrown.
//
bool throwAsyncJavaException( JNIEnv* env, const char* message ) {
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

// Re-throw a PgenException (that originated in underlying *pgen-lib* code) as a Java exception.
// Note that control RETURNS to the caller after the exception is thrown.
bool reThrowAsAsyncJavaException( JNIEnv* env, const PgenException& pgenException, const char* context) {
    static constexpr int kReservedMessageBufSize = 1024;
    static char reservedForExceptionMessage[kReservedMessageBufSize];
    snprintf(reservedForExceptionMessage, kReservedMessageBufSize, "%s / %s", pgenException.what(), context);
    return throwAsyncJavaException(env, pgenException.what());
}
