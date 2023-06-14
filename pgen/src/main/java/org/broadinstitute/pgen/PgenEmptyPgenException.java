package org.broadinstitute.pgen;

/**
 * Exception to signal that no variants have been written to the PGEN. Empty pgens are not valid.
 */
public class PgenEmptyPgenException extends RuntimeException {
    public PgenEmptyPgenException(final String message) {
        super(message);
    }
    
}
