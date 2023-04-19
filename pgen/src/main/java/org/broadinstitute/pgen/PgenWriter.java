/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package org.broadinstitute.pgen;

import htsjdk.io.HtsPath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import java.nio.BufferOverflowException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PgenWriter implements VariantContextWriter {

    public static final int NO_CALL_VALUE = -9;
    public static final int MAX_PLINK2_ALTERNATE_ALLELES = 255;

    public enum PgenWriteMode {
        PGEN_FILE_MODE_BACKWARD_SEEK(0),
        PGEN_FILE_MODE_WRITE_SEPARATE_INDEX(1),
        PGEN_FILE_MODE_WRITE_AND_COPY(2); 

        private final int mode;
        
        private PgenWriteMode(final int mode) {
            this.mode = mode;
        }

        public int value() { return this.mode; }
    };

    private long pgenContextHandle;
    private ByteBuffer alleleBuffer;

    private int maxAltAlleles = MAX_PLINK2_ALTERNATE_ALLELES;

    static {
        System.loadLibrary("pgen");
    }

    // doesn't preserve phasing
    public PgenWriter(final HtsPath file, final PgenWriteMode pgenWriteMode, final int maxAltAlleles, final long numberOfVariants, final int numberOfSamples) {
        if (maxAltAlleles > MAX_PLINK2_ALTERNATE_ALLELES) {
            throw new PgenJniException(
                String.format("Requested max alternate alleles of (%d) exceeds the supported pgen max of %d", maxAltAlleles, MAX_PLINK2_ALTERNATE_ALLELES));
        }
        this.maxAltAlleles = maxAltAlleles;

        pgenContextHandle = openPgen(file.getRawInputString(), pgenWriteMode.value(), numberOfVariants, numberOfSamples);
        alleleBuffer = createBuffer(numberOfSamples*2*4); //samples * ploidy * bytes in int32
        alleleBuffer.order(ByteOrder.LITTLE_ENDIAN);
    }

    @Override
    public void writeHeader(final VCFHeader header) {
       throw new UnsupportedOperationException("PGEN writer does not support independent header write.");
    }

    @Override
    public void setHeader(final VCFHeader header) {
        throw new UnsupportedOperationException("PGEN writer does not support independent setHeader");
   }

   @Override
    public void close() {
        closePgen(pgenContextHandle);
        pgenContextHandle = 0;
        destroyByteBuffer(alleleBuffer);
        alleleBuffer = null;
    }

    @Override
    public boolean checkError() {
        return false;
    }

    private static Map<Allele, Integer> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, Integer> alleleMap = new HashMap<>(vc.getAlleles().size() + 1);
        alleleMap.put(Allele.NO_CALL, NO_CALL_VALUE); // convenience for lookup
        final List<Allele> alleles = vc.getAlleles();
        for (int i = 0; i < alleles.size(); i++) {
            alleleMap.put(alleles.get(i), i);
        }

        return alleleMap;
    }

    @Override
    public void add(final VariantContext vc) {
        if (!vc.isBiallelic()) {
            throw new PgenJniException(
                String.format("Variant has %d alleles - multi-allelic variants are not yet implemented by the pgen-writer: %s",
                 vc.getNAlleles(),
                 vc.toStringWithoutGenotypes()));
        }

        //reset buffer
        alleleBuffer.clear();
        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);

        //System.out.println("\nVariant: " + vc.getContig() + "/" + vc.getStart());
        for (final Genotype g : vc.getGenotypes()) {
            if (g.getPloidy() != 2) {
                throw new PgenJniException(
                    "PGEN only supports diploid samples and we see one with ploidy = " + g.getPloidy()
                        + " at line " + vc.toStringDecodeGenotypes());
            }
            //System.out.println("  Genotype: " + g.getSampleName());
            for (final Allele allele : g.getAlleles()) {
                final Integer mapping = alleleMap.get(allele);
                //System.out.println("    Allele: " + allele.getBaseString() + ": " + mapping);
                try {
                    alleleBuffer.putInt(mapping);
                } catch (BufferOverflowException e){
                    throw new RuntimeException("Buffer overflow for: " + mapping +" for  Allele: " + allele.toString() + " from Genotype: " + g.toString() + " at buffer position: "+ alleleBuffer.position());
                }
            }
        }
        if (alleleBuffer.position() != alleleBuffer.limit()) {
            throw new IllegalStateException("Allele buffer is not completely filled, we have a problem. " +
                    "Position: " + alleleBuffer.position() + " Expected " + alleleBuffer.limit());
        }
        alleleBuffer.rewind();
        appendAlleles(pgenContextHandle, alleleBuffer);
    }

    private static native long openPgen(String file, int pgenWriteModeInt, long numberOfVariants, int numberOfSamples);
    // private static native void appendBiallelic(long pgenContextHandle, )
    private native void closePgen(long pgenContextHandle);
    private native void appendAlleles(long pgenContextHandle, ByteBuffer alleles);

    private static native ByteBuffer createBuffer(int length);
    private static native void destroyByteBuffer(ByteBuffer buffer);
}
