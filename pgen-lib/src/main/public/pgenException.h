//
//

#ifndef PGEN_LIB_PGENEXCEPTION_H
#define PGEN_LIB_PGENEXCEPTION_H
#include <exception>

//TODO: copy constructor (defaults to member-wise)
//TODO: free message - who owns it ?
//TODO: virtual method in subclass req'd?
//TODO: subclassing PgenException ?
//TODO: const syntax
//TODO: nullptr <-> NULL

// Exception class for passing results back to callers
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

    //    PgenException(const PgenException&) throw() {
    //    }

};

#endif //PGEN_LIB_PGENEXCEPTION_H
