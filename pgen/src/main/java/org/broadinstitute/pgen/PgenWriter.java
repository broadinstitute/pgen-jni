/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
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
import java.util.List;
import java.util.Map;

/**
 * An [HTSJDK](https://github.com/samtools/htsjdk) [VariantContextWriter]
 * (https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/writer/VariantContextWriter.java)
 * that writes [VariantContext](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/VariantContext.java)
 * objects to a file in the [plink2](https://www.cog-genomics.org/plink/2.0) in the plink2 PGEN format. See the plink2 PGEN
 * [spec](https://github.com/chrchang/plink-ng/tree/master/pgen_spec) for more information about the PGEN file format.
 * 
 * The writer implementation uses an underlying native component, which is built from a combination of local source files,
 * plus some source files taken directly from the plink2 (that are used by plink2 to build the pgen-lib target). Only
 * specific platforms are supported.
*/
public class PgenWriter implements VariantContextWriter {
    private static Log logger = Log.getInstance(PgenWriter.class);

    /**
     * A sentinel to signal that the number of variants is unknown. This value can be used as a variant count in the
     * pgen-lib APIs when the variant count is not known up front, as long as the corresponding file mode param is
     * not PGEN_FILE_MODE_BACKWARD_SEEK (PGEN_FILE_MODE_BACKWARD_SEEK requires an accurate variant count).
     */
    public static long VARIANT_COUNT_UNKNOWN = 0x7ffffffd; // plink2::kPglMaxVariantCt

    /**
     * Sentinel used by plink2 for missing data.
     */
    public static final int PLINK2_NO_CALL_VALUE = -9;
    /**
     * The maximum number of alternate alleles that plink2/pgen can handle (this is defined by plink2)
     */
    public static final int PLINK2_MAX_ALTERNATE_ALLELES = 254;  // plink2::kPglMaxAltAlleleCt

    public static String PGEN_EXTENSION = ".pgen";
    public static String PGEN_INDEX_EXTENSION = ".pgen.pgi";    
    public static String PVAR_EXTENSION = ".pvar";
    public static String PSAM_EXTENSION = ".psam";
 
    // If this java property is set/exists, the pgen native component will be loaded from the java.libary.path,
    // otherwise it is assumed to be included as a resource at the top level of a jar on the classpath.
    private static final String LOAD_PGEN_FROM_LIBRARY_PATH = "LOAD_PGEN_FROM_LIBRARY_PATH";

    /**
     * Enum for representing the plink2 pgen file write modes. See plink2::PgenWriteMode.
     */
    public enum PgenWriteMode {
        /**
         * requires backward seeks when writing
         */
        PGEN_FILE_MODE_BACKWARD_SEEK(0),
        /**
         * write a separate .pgi index file
         */
        PGEN_FILE_MODE_WRITE_SEPARATE_INDEX(1),
        /**
         * the final real .pgen is only created at the end, by writing the index and then appending the body of the first
         * temporary .pgen (which is then deleted).
         */
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

    // ******************** Native JNI methods  ********************
    private static native long openPgen(String file, int pgenWriteModeInt, long numberOfVariants, int numberOfSamples, int maxAltAlleles);
    private static native boolean closePgen(long pgenContextHandle, long numDroppedVariants);
    private static native long getPgenVariantCount(long pgenContextHandle);
    private static native boolean appendAlleles(long pgenContextHandle, ByteBuffer alleles, int alleleCount);
    private static native ByteBuffer createBuffer(int length);
    private static native boolean destroyByteBuffer(ByteBuffer buffer);
   // ******************** End Native JNI methods  ********************
 
    static {
        if (System.getProperty(LOAD_PGEN_FROM_LIBRARY_PATH) != null) {
            // for local testing within the IDE
            System.loadLibrary("pgen");
        } else {
            // otherwise, load it from a jar file on the classpath
            NativeLibraryUtils.loadLibraryFromClasspath(
                NativeLibraryUtils.runningOnMac() ?
                    "/libpgen.dylib" :
                    "/libpgen.so"
            );
        }
    }
 
    /**
     * Create a PGEN writer. The writer creats a PGEN file set (.pgen and .pvar/.psam files). Depending on the file
     * mode used, may also create a .pgen.pgi file.
     * 
     * Supports only diploid sites with fewer than {@code PLINK2_MAX_ALTERNATE_ALLELES} (this value can be
     * furhter limited using {@code maxAltAlleles}. Sites that exceed this number are silently dropped.
     * 
     * NOTE: Currently this writer doesn't preserve genotype phasing.
     */
    public PgenWriter(
        final HtsPath pgenFile,
        final VCFHeader vcfHeader,
        final PgenWriteMode pgenWriteMode,
        final long numberOfVariants,
        final int maxAltAlleles) {

        if (!pgenFile.hasExtension(PGEN_EXTENSION)) {
            throw new PgenException(
                String.format("Invalid pgen file name: %s. pgen files must use the .pgen extension", pgenFile.getRawInputString()));
        }
        if (!pgenFile.getScheme().equals("file")) {
            throw new PgenException(String.format("Invalid pgen file name: %s. pgen files must be local files", pgenFile));
        }
        if (maxAltAlleles > PLINK2_MAX_ALTERNATE_ALLELES) {
            throw new PgenException(
                String.format("Requested max alternate alleles of (%d) exceeds the supported pgen max of (%d)",
                    maxAltAlleles,
                    PLINK2_MAX_ALTERNATE_ALLELES));
        }
        this.maxAltAlleles = maxAltAlleles;
        this.expectedVariantCount = numberOfVariants;

        pgenContextHandle = openPgen(pgenFile.getRawInputString(), pgenWriteMode.value(), numberOfVariants, vcfHeader.getNGenotypeSamples(), maxAltAlleles);
        if (pgenContextHandle == 0) {
            //openPgen threw an async Java exception
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

        // closePgen returns false if it had to throw an async Java exception, so test for that, and if it failed,
        // don't do anything else that might throw.
        //
        // Tell the writer how many variants we dropped (due to exceeding the # of alternate alleles) so it
        // doesn't throw if the number written doesn't match the number expected (which is provided when the
        // writer is opened)
       if (closePgen(pgenContextHandle, droppedVariantCount)) {
            pgenContextHandle = 0;
            //destroyByteBuffer might return false if for some reason it has to throw an async Java exception, but
            // we don't need to test for that here since we're only nulling out a variable on return
            destroyByteBuffer(alleleBuffer);
            alleleBuffer = null;    
       }
    }

    @Override
    public boolean checkError() {
        return false;
    }

    @Override
    public void add(final VariantContext vc) {
        if (vc.getNAlleles() > maxAltAlleles) {
            droppedVariantCount++;
            return;
        }
        
        //reset buffer
        alleleBuffer.clear();
        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);
     
        for (final Genotype g : vc.getGenotypes()) {
            if (g.getPloidy() != 2) {
                throw new PgenException(
                    "PGEN only supports diploid samples and we see one with ploidy = " + g.getPloidy()
                        + " at line " + vc.toStringWithoutGenotypes());
            }
           for (final Allele allele : g.getAlleles()) {
                final Integer mapping = alleleMap.get(allele);
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
        final boolean appendRet = appendAlleles(pgenContextHandle, alleleBuffer, alleleMap.size() - 1);
        
        if (appendRet) {
            // only add to the pvar if appendAlleles succeeded
            pVarWriter.add(vc);
        }
    }

   /**
     * @return the number of variants dropped because they exceeded the max alternate allele count
     */
    public long getDroppedVariantCount() { return droppedVariantCount; }

     /**
     * @return the number of variants actually written to the pgen
     * 
     * Delegates to the pgen-lib code to get the actual number recorded by the pgen library code.
     */
    public long getWrittenVariantCount() { return getPgenVariantCount(pgenContextHandle); }

    /**
     * given a Path, return the absolute path of the file, without the trailing extension
     */
    // Visible for testing
    public static String getAbsoluteFileNameWithoutExtension(final Path targetPath, final String extension) {
        final String targetAbsolutePath = targetPath.toAbsolutePath().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf(extension));
    }
    
    /**
     * Create a .pvar companion file for {@code pgenFile}.
     */
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

    /**
     * Create and writes a .psam companion file for {@code pgenFile}.
     */
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
