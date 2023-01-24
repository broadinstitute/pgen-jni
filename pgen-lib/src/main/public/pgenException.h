//
//

#ifndef PGEN_LIB_PGENEXCEPTION_H
#define PGEN_LIB_PGENEXCEPTION_H
#include <exception>
#include "plink2_base.h"

namespace pgenlib {
    // reserve some static memory to be used for consructing error messages for exceptions
    static constexpr int kReservedMessageBufSize = 1024;
    static char reservedForExceptionMessage[kReservedMessageBufSize];

    // Exception class for passing results back to pgenlib callers
    class PgenException : public std::exception {
        private:
            const char *message;

        public:
            PgenException(const char* message) {
                this->message = message;
            }

            virtual const char* what() const throw() {
                return message;
            }

            PgenException(const PgenException&) throw() = default;
    };
}

#endif //PGEN_LIB_PGENEXCEPTION_H
