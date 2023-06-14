//
//

#ifndef PGEN_LIB_PGENVARIANTCOUNTEXCEPTION_H
#define PGEN_LIB_PGENVARIANTCOUNTEXCEPTION_H
#include <exception>
#include "plink2_base.h"

namespace pgenlib {
    // Exception class for the specific case where too few variants have been written when closePgen is called
    class PgenVariantCountException : public std::exception {
    private:
        const char *message;

    public:
        PgenVariantCountException(const char* message) {
            // make a copy in our reserved memory, since the caller is probably about to throw...
            this->message = strncpy(reservedForExceptionMessage, message, kReservedMessageBufSize);
        }

        virtual const char* what() const throw() {
            return message;
        }

        PgenVariantCountException(const PgenVariantCountException&) throw() = default;
    };
}

#endif //PGEN_LIB_PGENVARIANTCOUNTEXCEPTION_H
