#ifndef __PGEN_JNI_UTILS__
#define __PGEN_JNI_UTILS__

#include "pgenlib_write.h"
#include <jni.h>

//Check and throw a java exception if we returned an error.
void checkPglErr(JNIEnv* env, const plink2::PglErr result, const char* message){
    if (result){ // PglErr evaluates to true if it was an error
        jclass myExceptionClass = env->FindClass("org/broadinstitute/pgen/PgenJniException");
        jstring myErrorJString = env->NewStringUTF(message);
        jmethodID ctorMethod = env->GetMethodID(myExceptionClass, "<init>", "(Ljava/lang/String;I)V");
        jthrowable myExceptionObject = (jthrowable)env->NewObject(myExceptionClass, ctorMethod, myErrorJString, (jint)(uint32_t)result);
        env->Throw(myExceptionObject);
    }
}
#endif // __PGEN_JNI_UTILS__