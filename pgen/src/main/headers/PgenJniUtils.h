#ifndef __PGEN_JNI_UTILS_H__
#define __PGEN_JNI_UTILS_H__

#include "pgenException.h"
#include <jni.h>

using namespace pgenlib;

// Throw a Java exception (PgenJniException) with the given error message. Note that control RETURNS
// to the caller after the exception is thrown.
//
bool throwAsyncJavaException( JNIEnv* env, const char* message, const char *javaExceptionClassName );

// Re-throw a PgenException (that originated in underlying *PGENLIB* code) as a Java exception.
//  Note that control RETURNS to the caller after the exception is thrown.
bool reThrowAsAsyncJavaException( JNIEnv* env, const PgenException& pgenException, const char* context);

#endif // __PGEN_JNI_UTILS_H__