/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package org.broadinstitute.pgen;

import htsjdk.io.HtsPath;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFHeader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.BufferOverflowException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.HashMap;
//import java.util.HashSet;
import java.util.List;
import java.util.Map;
//import java.util.Set;

public class PgenWriter implements VariantContextWriter {
    private static Log logger = Log.getInstance(PgenWriter.class);

    public static final int PLINK2_NO_CALL_VALUE = -9;
    public static final int PLINK2_MAX_ALTERNATE_ALLELES = 254;  // plink2::kPglMaxAltAlleleCt

    public static String PGEN_EXTENSION = ".pgen";
    public static String PGEN_INDEX_EXTENSION = ".pgen.pgi";    
    public static String PVAR_EXTENSION = ".pvar";
    public static String PSAM_EXTENSION = ".psam";
 
    // enum for the plink2 pgen write modes
    public enum PgenWriteMode {
        PGEN_FILE_MODE_BACKWARD_SEEK(0),
        PGEN_FILE_MODE_WRITE_SEPARATE_INDEX(1),
        PGEN_FILE_MODE_WRITE_AND_COPY(2); 

        private final int mode;
        private PgenWriteMode(final int mode) { this.mode = mode; }
        public int value() { return this.mode; }
    };

    private HtsPath pVarFile = null;
    private HtsPath pSamFile = null;
    private final int maxAltAlleles;
    private long pgenContextHandle;
    private ByteBuffer alleleBuffer;
    private VariantContextWriter pVarWriter;
    private long expectedVariantCount = 0L;
    private long droppedVariantCount = 0L;
    // private long multiallelic_ct = 0;
    // private long nonSNP_ct = 0;
    // private long mnp_ct = 0;

    // ******************** Native JNI methods  ********************
    private static native long openPgen(String file, int pgenWriteModeInt, long numberOfVariants, int numberOfSamples, int maxAltAlleles);
    private static native boolean closePgen(long pgenContextHandle, long numDroppedVariants);
    /**
     * @return the number of variants actually written to the pgen
     */
    private static native long getPgenVariantCount(long pgenContextHandle);
    private static native boolean appendAlleles(long pgenContextHandle, ByteBuffer alleles, int alleleCount);
    private static native ByteBuffer createBuffer(int length);
    private static native boolean destroyByteBuffer(ByteBuffer buffer);
   // ******************** End Native JNI methods  ********************
 
    static {
        System.loadLibrary("pgen");
    }

    // doesn't preserve phasing
    public PgenWriter(
        final HtsPath pgenFile,
        final VCFHeader vcfHeader,
        final PgenWriteMode pgenWriteMode,
        final long numberOfVariants,
        final int maxAltAlleles) {

        if (!pgenFile.hasExtension(PGEN_EXTENSION)) {
            throw new PgenJniException(
                String.format("Invalid pgen file name: %s. pgen files must use the .pgen extension", pgenFile.getRawInputString()));
        }
        if (!pgenFile.getScheme().equals("file")) {
            throw new PgenJniException(String.format("Invalid pgen file name: %s. pgen files must be local files", pgenFile));
        }
        if (maxAltAlleles > PLINK2_MAX_ALTERNATE_ALLELES) {
            throw new PgenJniException(
                String.format("Requested max alternate alleles of (%d) exceeds the supported pgen max of %d",
                    maxAltAlleles,
                    PLINK2_MAX_ALTERNATE_ALLELES));
        }
        this.maxAltAlleles = maxAltAlleles;
        this.expectedVariantCount = numberOfVariants;

        pgenContextHandle = openPgen(pgenFile.getRawInputString(), pgenWriteMode.value(), numberOfVariants, vcfHeader.getNGenotypeSamples(), maxAltAlleles);
        if (pgenContextHandle == 0) {
            //openPgen wthrew an async Java exception
            return;
        }
        alleleBuffer = createBuffer(vcfHeader.getNGenotypeSamples() * 2 * 4); //samples * ploidy * bytes in int32_t (sizeof AlleleCode)
        if (alleleBuffer == null) {
            //createBuffer threw an async Java exception
            return;
        }
        alleleBuffer.order(ByteOrder.LITTLE_ENDIAN);
 
        // create the .pvar, and write the entire psam
        pVarFile = createPVAR(pgenFile, vcfHeader);
        pSamFile = writePSAM(pgenFile, vcfHeader);
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
        //System.out.println(String.format("Multiallelic: %d NonSNP: %d MNP: %d", multiallelic_ct, nonSNP_ct, mnp_ct));
   
        pVarWriter.close();
        pVarWriter = null;

        // tell the writer how many variants we dropped (due to exceeding the # of alternate alleles) so it
        // doesn't throw if the number written doesn't match the number expected (which is provided when the
        // writer is opened)
        closePgen(pgenContextHandle, droppedVariantCount);
        pgenContextHandle = 0;
        destroyByteBuffer(alleleBuffer);
        alleleBuffer = null;
    }

    @Override
    public boolean checkError() {
        return false;
    }

    @Override
    public void add(final VariantContext vc) {
        // if (!vc.isBiallelic()) {
        //     multiallelic_ct++;
        // }
        // if (!vc.isSNP()) {
        //     nonSNP_ct++;
        // }
        // if (vc.isMNP()) {
        //     mnp_ct++;
        // }

        if (vc.getNAlleles() > maxAltAlleles) {
            droppedVariantCount++;
            return;
        }
        
        //reset buffer
        alleleBuffer.clear();
        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);

//        final Set<Integer> observedAlleles = new HashSet<>();
        for (final Genotype g : vc.getGenotypes()) {
            if (g.getPloidy() != 2) {
                throw new PgenJniException(
                    "PGEN only supports diploid samples and we see one with ploidy = " + g.getPloidy()
                        + " at line " + vc.toStringWithoutGenotypes());
            }
           for (final Allele allele : g.getAlleles()) {
                final Integer mapping = alleleMap.get(allele);
                // if (mapping != -9) {
                //     // don't count missing dat as an observed allele
                //     observedAlleles.add(mapping);
                // }
                try {
                    alleleBuffer.putInt(mapping);
                } catch (final BufferOverflowException e) {
                    throw new RuntimeException(
                        String.format(
                            "Buffer overflow at position: %d mapping: %d for variant: %s, genotype: %s allele: %s",
                            alleleBuffer.position(),
                            mapping,
                            vc.toStringWithoutGenotypes(),
                            g.toString(),
                            allele.toString()),
                        e);
                }
            }
        }
        if (alleleBuffer.position() != alleleBuffer.limit()) {
            throw new IllegalStateException("Allele buffer is not completely filled, we have a problem. " +
                    "Position: " + alleleBuffer.position() + " Expected " + alleleBuffer.limit());
        }
        alleleBuffer.rewind();
        // -1 to account for the synthetic -9 value that is always in in the alleleMap
        int alleleMapSize_actual = alleleMap.size() - 1;
        // if (observedAlleles.size() != alleleMapSize_actual) {
        //     System.out.println(String.format(
        //         "Alleles not observed (%d/%d) pos: %d",
        //          observedAlleles.size(),
        //         alleleMapSize_actual,
        //          vc.getStart()));
        // }
        final boolean appendRet = appendAlleles(pgenContextHandle, alleleBuffer, alleleMap.size() - 1);
        if (appendRet) {
            // only add to the pvar if appendAlleles succeeded
            pVarWriter.add(vc);
        }
    }

   /**
     * @return the number of variants dropped due to exceeding the max alternate allele count
     */
    public long getDroppedVariantCount() { return droppedVariantCount; }

     /**
     * @return the number of variants actually written to the pgen
     * 
     * Delegate to the pgen-lib code to get the actual number recorded by the pgen library code.
     */
    public long getWrittenVariantCount() { return getPgenVariantCount(pgenContextHandle); }

    // given a Path, return the absolute path of the file, without the trailing extension
    public static String getAbsoluteFileNameWithoutExtension(final Path targetPath, final String extension) {
        final String targetAbsolutePath = targetPath.toAbsolutePath().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf(extension));
    }
    
    // create the .pvar
    private HtsPath createPVAR(final HtsPath pgenFile, final VCFHeader vcfHeader) {
        final String pgenFilePrefix = getAbsoluteFileNameWithoutExtension(pgenFile.toPath(), PGEN_EXTENSION);
        final HtsPath pVarFile = new HtsPath(pgenFile.toPath()
            .resolveSibling(pgenFilePrefix + PgenWriter.PVAR_EXTENSION)
            .toAbsolutePath().toString());
        pVarWriter = new VariantContextWriterBuilder()
            .clearOptions()
            .setOptions(EnumSet.of(Options.DO_NOT_WRITE_GENOTYPES, Options.ALLOW_MISSING_FIELDS_IN_HEADER))
            .setOutputPath(pVarFile.toPath())
            .setOutputFileType(OutputType.VCF) // plink2 expects the .pvar to have a .pvar extension
            .build();
        pVarWriter.writeHeader(vcfHeader);
        return pVarFile;
    }

    // write the entire psam out
    private HtsPath writePSAM(final HtsPath pgenFile, final VCFHeader vcfHeader) {
        final String PSAM_HEADER_LINE = "#IID\tSEX\n";
        final String PSAM_DETAIL_LINE = "\tN/A\n";

        // create, write, and close the .psam up front, so we don't have to retain the header until the end
        final String pgenFilePrefix = getAbsoluteFileNameWithoutExtension(pgenFile.toPath(), PGEN_EXTENSION);
        final HtsPath pSamFile = new HtsPath(pgenFile.toPath()
                        .resolveSibling(pgenFilePrefix + PgenWriter.PSAM_EXTENSION)
                        .toAbsolutePath().toString());
        try (final BufferedWriter psamWriter = Files.newBufferedWriter(pSamFile.toPath())) {
            psamWriter.append(PSAM_HEADER_LINE);
            // Sample name order matters here. I don't see any spec for the psam or pvar, but if you use plink2 to crecrate a VCF
            // from a pgen file set, it appears to use the order of the samples in the .psam as the basis for linking the genotypes
            // in the pgen back to the VCF. So if we don't preserve the order in the .psam, the genotypes in the roundtripped VCF
            // won't match the original VCF (and will be incorrect).
            for (final String sampleName : vcfHeader.getGenotypeSamples()) {
                psamWriter.write(sampleName);
                psamWriter.write(PSAM_DETAIL_LINE);
            }
        } catch (final IOException e) {
            throw new RuntimeIOException(String.format("Error writing the .psam file %s", pSamFile.getRawInputString()), e);
        }
        return pSamFile;
    }

    private static Map<Allele, Integer> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, Integer> alleleMap = new HashMap<>(vc.getAlleles().size() + 1);
        alleleMap.put(Allele.NO_CALL, PLINK2_NO_CALL_VALUE); // convenience for lookup
        final List<Allele> alleles = vc.getAlleles();
        for (int i = 0; i < alleles.size(); i++) {
            alleleMap.put(alleles.get(i), i);
        }

        return alleleMap;
    }

}
