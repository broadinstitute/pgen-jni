package org.broadinstitute.pgen;

public class PgenJniException extends RuntimeException {
    final String message;

    public PgenJniException(final String message){
        super(message);
        this.message = message;
    }

}
