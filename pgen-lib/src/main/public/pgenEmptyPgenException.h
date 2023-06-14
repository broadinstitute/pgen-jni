//
//

#ifndef PGEN_LIB_PGENEMPTYPGENEXCEPTION_H
#define PGEN_LIB_PGENEMPTYPGENEXCEPTION_H
#include <exception>
#include "pgenException.h"

namespace pgenlib {

    // Exception class for the specific case where No Ovariants have been written when closePgen is called
    class PgenEmptyPgenException : public std::exception {
    private:
        const char *message;

    public:
        PgenEmptyPgenException(const char *message) {
            // make a copy in our reserved memory, since the caller is probably about to throw...
            this->message = strncpy(reservedForExceptionMessage, message, kReservedMessageBufSize);
        }

        virtual const char *what() const throw() {
            return message;
        }

        PgenEmptyPgenException(const PgenEmptyPgenException &) throw() = default;
    };
}
#endif //PGEN_LIB_PGENEMPTYPGENEXCEPTION_H
