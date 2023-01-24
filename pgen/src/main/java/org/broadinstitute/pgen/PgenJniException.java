package org.broadinstitute.pgen;

public class PgenJniException extends RuntimeException {
    final String message;

    public PgenJniException(final String message){
        super(message);
        this.message = message;
    }

    // public PgenJniException(String message, int errorCode){
    //     this(message + " caused by: " + PglErr.getCorrespondingError(errorCode) + "(" + errorCode +")");
    // }

    // public enum PglErr {
    //     Success(0),
    //     Skipped(1),
    //     Nomem(2),
    //     OpenFail(3),
    //     ReadFail(4),
    //     WriteFail(5),
    //     // MalformedInput should be returned on low-level file format violations,
    //     // while InconsistentInput should be returned for higher-level logical
    //     // problems like mismatched files (generally solvable by fixing the command
    //     // line), and DegenerateData for properly-formatted-and-matched files that
    //     // yields degenerate computational results due to e.g. divide by zero or
    //     // insufficient rank.
    //     MalformedInput(6),
    //     InconsistentInput(7),
    //     InvalidCmdline(8),
    //     ThreadCreateFail(9),
    //     NetworkFail(10),
    //     VarRecordTooLarge(11),
    //     UnsupportedInstructions(12),
    //     DegenerateData(13),
    //     DecompressFail(14), // also distinguish this from MalformedInput
    //     RewindFail(15),
    //     GpuFail(16),
    //     SampleMajorBed(32),
    //     NomemCustomMsg(59),
    //     InternalError(60),
    //     WarningErrcode(61),
    //     ImproperFunctionCall(62),
    //     NotYetSupported(63),
    //     // These are only for internal use.  If any of these reach the top level
    //     // instead of being handled or converted to another error code, that's a bug,
    //     // and plink2 prints a message to that effect.
    //     Help(125),
    //     LongLine(126),
    //     Eof(127),

    //     //this is what we get if we get back an error code I haven't captured here.
    //     ErrorCodeNotFound(Integer.MIN_VALUE);
    //     private final int errorCode;

    //     PglErr(int errorCode){
    //         this.errorCode = errorCode;
    //     }

    //     private static final Map<Integer, PglErr> errorCodeToValue = new LinkedHashMap<>();

    //     static {
    //         for( PglErr err: PglErr.values()) {
    //             errorCodeToValue.put(err.errorCode, err);
    //         }
    //     }

    //     public static PglErr getCorrespondingError(int errorCode){
    //         return errorCodeToValue.getOrDefault(errorCode, ErrorCodeNotFound);
    //     }

    //     public int getErrorCode(){
    //         return errorCode;
    //     }
    // }
}
