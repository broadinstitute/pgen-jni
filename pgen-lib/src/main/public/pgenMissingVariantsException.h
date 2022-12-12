//
//

#ifndef PGEN_LIB_PGENMISSINGVARIANTSEXCEPTION_H
#define PGEN_LIB_PGENMISSINGVARIANTSEXCEPTION_H
#include <exception>
#include "pgenException.h"

namespace pgenlib {
    // Exception class for the specific case where too few variants have been written when closePgen is called
    class PgenMissingVariantsException : public std::exception {
    private:
        const char *message;

    public:
        PgenMissingVariantsException(const char* message) {
            // make a copy in our reserved memory, since the caller is probably about to throw...
            this->message = strncpy(reservedForExceptionMessage, message, kReservedMessageBufSize);
        }

        virtual const char* what() const throw() {
            return message;
        }

        PgenMissingVariantsException(const PgenMissingVariantsException&) throw() = default;
    };
}

#endif //PGEN_LIB_PGENMISSINGVARIANTSEXCEPTION_H
