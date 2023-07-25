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
 * plus some source files taken directly from the plink2 (that are used by plink2 to build the pgenlib target). Only
 * some platforms are supported (macos and linux intel platforms).
*/
public class PgenWriter implements VariantContextWriter {
    private static Log logger = Log.getInstance(PgenWriter.class);

    /**
     * variant count is not known up front, as long as the corresponding file mode param used is either {@code #PGEN_FILE_MODE_WRITE_SEPARATE_INDEX}
     * or {@code #PGEN_FILE_MODE_WRITE_AND_COPY} (the write mode must not be {@link #PGEN_FILE_MODE_BACKWARD_SEEK}, which requires an accurate
     * upfront variant count).
     */
    public static long VARIANT_COUNT_UNKNOWN = 0x7ffffffd; // plink2::kPglMaxVariantCt

    /**
     * Sentinel used by plink2 for missing data.
     */
    public static final int PLINK2_NO_CALL_VALUE = -9;
    /**
     * The maximum number of alternate alleles that plink2/PGEN can handle (this is determined/defined by plink2)
     */
    public static final int PLINK2_MAX_ALTERNATE_ALLELES = 254;  // plink2::kPglMaxAltAlleleCt

    public static String PGEN_EXTENSION = ".pgen";
    public static String PGEN_INDEX_EXTENSION = ".pgen.pgi";    
    public static String PVAR_EXTENSION = ".pvar";
    public static String PSAM_EXTENSION = ".psam";
 
    /**
     * Enum for representing the plink2 PGEN file write modes. See plink2::PgenWriteMode.
     */
    public enum PgenWriteMode {
        /**
         * write the pgen in single pass; requires backward seeks when writing; can only be used when an accurate upfront variant count
         * is provided to the {@link PgenWriter}
         */
        PGEN_FILE_MODE_BACKWARD_SEEK(0),        // plink2::PgenWriteMode::kPgenWriteBackwardSeek
        /**
         * write a separate .pgi index file
         */
        PGEN_FILE_MODE_WRITE_SEPARATE_INDEX(1), // plink2::PgenWriteMode::kPgenWriteSeparateIndex
        /**
         * the final real .pgen is only created at the end, by writing the index and then appending the body of the first
         * temporary .pgen (which is then deleted).
         */
        PGEN_FILE_MODE_WRITE_AND_COPY(2);       // plink2::PgenWriteMode::kPgenWriteAndCopy
 
        private final int mode;
        private PgenWriteMode(final int mode) { this.mode = mode; }
        public int value() { return this.mode; }
    };

    /**
     * Enum for supported write mode flags.
     */
    public enum PgenWriteFlag {
        // This enum, and the corresponding enum values must be kept in sync with the corresponding constants
        // in pgenlib::PgenWriteFlags.
        PRESERVE_PHASING(0x1),  // pgenlib::kWriteFlagPreservePhasing
        MULTI_ALLELIC(0x2);     // pgenlib::kWriteFlagMultiAllelic

        private final int flag;
        private PgenWriteFlag(final int flag) { this.flag = flag; }
        public int value() { return this.flag; }

        /**
         * Convert an EnumSet<PgenWriteFlag> into the corresponding pgenlib bitwise/integer flags.
         */
        private static int toIntFlags(final EnumSet<PgenWriteFlag> flagsSet) {
            return (flagsSet.contains(PRESERVE_PHASING) ? PRESERVE_PHASING.value() : 0) |
                   (flagsSet.contains(MULTI_ALLELIC) ? MULTI_ALLELIC.value() : 0);
        }
    }

    private static byte PHASED_CODE = (byte) 1;
    private static byte UNPHASED_CODE = (byte) 0;

    private List<String> sampleNames = null;
    private HtsPath pVarFile = null;
    private HtsPath pSamFile = null;
    private final int maxAltAlleles;
    private long pgenContextHandle;
    private ByteBuffer alleleBuffer;
    private ByteBuffer phasingBuffer;
    private VariantContextWriter pVarWriter;
    private long expectedVariantCount = 0L;
    private long droppedVariantCount = 0L;

    // ******************** Native JNI methods  ********************
    private static native long openPgen(String file, int pgenWriteModeInt, int writeFlags, long numberOfVariants, int numberOfSamples, int maxAltAlleles);
    private static native boolean closePgen(long pgenContextHandle, long numDroppedVariants);
    private static native long getPgenVariantCount(long pgenContextHandle);
    private static native boolean appendAlleles(long pgenContextHandle, ByteBuffer alleles, ByteBuffer phasing, int alleleCount);
    private static native ByteBuffer createBuffer(int length);
    private static native boolean destroyByteBuffer(ByteBuffer buffer);
   // ******************** End Native JNI methods  ********************
 
    // If this java property is set/exists, the native PGEN writer component will be loaded from java.libary.path,
    // otherwise it is assumed to be included as a resource at the top level of a jar on the classpath.
    private static final String LOAD_PGEN_FROM_LIBRARY_PATH = "LOAD_PGEN_FROM_LIBRARY_PATH";
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
     * Create a PGEN writer for writing VariantContext objects to a Plink2 PGEN file. The writer creates a set of related PGEN
     * files (.pgen/.pvar/.psam files). Depending on the {@link #PgenWriteMode} specified, may also create a .pgen.pgi file,
     * or temporary/intermediate files.
     * 
     * Supports only diploid sites with fewer than {@link #PLINK2_MAX_ALTERNATE_ALLELES} (this value can be further limited using
     * {@code maxAltAlleles}). Sites that exceed the maximum number of alternate alleles are silently dropped.
     * 
     * @param pgenFileName the name of the PGEN file to be created (must have a ".pgen" suffix)
     * @param vcfHeader a valid VCF header for to use to create the pgen/pvar/psam files
     * @param pgenWriteMode the PGEN write mode to use (see {@link PgenWriteMode})
     * @param writeFlags the write flags to use - see {@link #PgenWriteFlag}. If phase information is present for the source genotypes, include
     * the {@link PgenWriteFlag#PRESERVE_PHASING} flag. If multi allelic variants are present, include the {@link PgenWriteFlag#MULTI_ALLELIC} flag.
     * @param maxAltAlleles the maximum number of alternate alleles to consider; site with more alterate alleles than this value will be
     * silently dropped
     */
    public PgenWriter(
        final HtsPath pgenFileName,
        final VCFHeader vcfHeader,
        final PgenWriteMode pgenWriteMode,
        final EnumSet<PgenWriteFlag> writeFlags,
        final int maxAltAlleles) {
            this(pgenFileName, vcfHeader, pgenWriteMode, writeFlags, VARIANT_COUNT_UNKNOWN, PLINK2_MAX_ALTERNATE_ALLELES);
    }

    /**
     * Create a PGEN writer for writing VariantContexts to a Plink2 PGEN file. The writer creates a set of related PGEN files (.pgen and
     * .pvar/.psam files). Depending on the {@link #PgenWriteMode} specified, may also create a .pgen.pgi file.
     * 
     * Supports only diploid sites with fewer than {@link #PLINK2_MAX_ALTERNATE_ALLELES} (this value can be further limited using
     * {@code maxAltAlleles}. Sites that exceed the maximum number of alternate alleles are silently dropped.
     * 
     * @param pgenFileName the name of the PGEN file to be created (must end in .pgen)
     * @param vcfHeader a valid VCF header to use to create the PGEN
     * @param pgenWriteMode the pgen write mode to use (see {@link PgenWriteMode#})
     * @param writeFlags the write flags to use - see {@link #PgenWriteFlag}. If phase information is present for the source genotypes, include
     * the {@link PgenWriteFlag#PRESERVE_PHASING} flag. If multi allelic variants are present, include the {@link PgenWriteFlag#MULTI_ALLELIC} flag.
     * @param numberOfVariants the number of variants to be written. if the number is unknown, use the sentinel value {@link PgenWriter#VARIANT_COUNT_UNKNOWN},
     * but doing so precludes the use of the PGEN write mode {@link PgenWriteMode#PGEN_FILE_MODE_BACKWARD_SEEK}
     * @param maxAltAlleles the maximum number of alternate alleles to consider; site with more alterate alleles than this value will be
     * silently dropped
     */
    public PgenWriter(
        final HtsPath pgenFileName,
        final VCFHeader vcfHeader,
        final PgenWriteMode pgenWriteMode,
        final EnumSet<PgenWriteFlag> writeFlags,
        final long numberOfVariants,
        final int maxAltAlleles) {

        if (!pgenFileName.hasExtension(PGEN_EXTENSION)) {
            throw new PgenException(
                String.format("Invalid PGEN file name: %s. PGEN files must use the .pgen extension", pgenFileName.getRawInputString()));
        }
        if (!pgenFileName.getScheme().equals("file")) {
            throw new PgenException(String.format("Invalid PGEN file name: %s. PGEN files must be local files", pgenFileName));
        }
        if (maxAltAlleles > PLINK2_MAX_ALTERNATE_ALLELES) {
            throw new PgenException(
                String.format("Requested max alternate alleles of (%d) exceeds the supported PGEN max of (%d)",
                    maxAltAlleles,
                    PLINK2_MAX_ALTERNATE_ALLELES));
        }
        this.maxAltAlleles = maxAltAlleles;
        this.expectedVariantCount = numberOfVariants;
        this.sampleNames = vcfHeader.getGenotypeSamples();

       pgenContextHandle = openPgen(
            pgenFileName.getRawInputString(),
            pgenWriteMode.value(),
            PgenWriteFlag.toIntFlags(writeFlags),
            numberOfVariants,
            vcfHeader.getNGenotypeSamples(),
            maxAltAlleles);
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
 
        phasingBuffer = createBuffer(vcfHeader.getNGenotypeSamples());
        if (phasingBuffer == null) {
            //createBuffer threw an async Java exception
            return;
        }
        phasingBuffer.order(ByteOrder.LITTLE_ENDIAN);

        // create the .pvar, and write the entire psam
        pVarFile = createPVAR(pgenFileName, vcfHeader);
        pSamFile = writePSAM(pgenFileName, vcfHeader);
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
            destroyByteBuffer(phasingBuffer);
            phasingBuffer = null;    
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
        
        alleleBuffer.clear();
        phasingBuffer.clear();
        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);
    
        // because there may be missing genotyopes, it is significantly cleaner code-wise to iterate through the 
        // sample names than through the genotypes, but this code does not recognize if there is a genotype in the VC
        // that has a sample name that is not in the header
        for (final String sampleName : sampleNames) {
            final Genotype g = vc.getGenotype(sampleName);
            if (g != null) {
                if (g.getPloidy() != 2) {
                    throw new PgenException(
                        String.format("PGEN only supports diploid samples but a sample with ploidy = %d was found at variant %s",
                            g.getPloidy(),
                            vc.toStringWithoutGenotypes()));
                } else if (g.getAlleles().size() != 2) {
                    throw new IllegalArgumentException(String.format("Bad allele count in genotype %d", g.getAlleles().size()));
                }
                for (final Allele allele : g.getAlleles()) {
                    //TODO: should this check for/handle the case where the allele is not in the allele map for the VC ?
                    updateAlleleBuffer(vc, g, allele, alleleMap.get(allele));
                }
                updatePhasingBuffer(vc, g, g.isPhased() ? PHASED_CODE : UNPHASED_CODE);
            } else {
                // synthesize no-call values for any missing genotypes
                updateAlleleBuffer(vc, null, null, PLINK2_NO_CALL_VALUE);
                updateAlleleBuffer(vc, null, null, PLINK2_NO_CALL_VALUE);
                updatePhasingBuffer(vc, null, UNPHASED_CODE);
           }
        }

         if (alleleBuffer.position() != alleleBuffer.limit()) {
            throw new IllegalStateException(
                String.format("Allele buffer is not completely filled, position is %d but expected %d.",
                    alleleBuffer.position(),
                    alleleBuffer.limit()));
        } else if (phasingBuffer.position() != phasingBuffer.limit()) {
            throw new IllegalStateException(
                String.format("Phase buffer is not completely filled, position is %d but expected %d.",
                    phasingBuffer.position(),
                    phasingBuffer.limit()));
        }

        alleleBuffer.rewind();
        phasingBuffer.rewind();
        final boolean appendRet = appendAlleles(pgenContextHandle, alleleBuffer, phasingBuffer, alleleMap.size() - 1);
        if (appendRet) { // only add to the pvar if appendAlleles succeeded
            pVarWriter.add(vc);
        }
    }

   /**
     * @return the number of variants dropped because they exceeded the max alternate allele count
     */
    public long getDroppedVariantCount() { return droppedVariantCount; }

     /**
     * @return the number of variants actually written to the PGEN
     * 
     * Delegates to the pgenlib code to get the actual number recorded by the pgen library code.
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
     * Creates writes a .psam companion file for {@code pgenFile}.
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
            // Sample name order matters here. If you use plink2 to create a VCF from a PGEN file set, it appears to use the order of the samples in
            // the .psam as the basis for linking the genotypes in the PGEN back to the VCF samples. So if we don't preserve the order in the .psam,
            // the genotypes in the roundtripped VCF won't match the original VCF, and will be incorrect.
            for (final String sampleName : vcfHeader.getGenotypeSamples()) {
                psamWriter.write(sampleName);
                psamWriter.write(PSAM_DETAIL_LINE);
            }
        } catch (final IOException e) {
            throw new RuntimeIOException(String.format("Error writing the .psam file %s", pSamFile.getRawInputString()), e);
        }
        return pSamFile;
    }

    private void updateAlleleBuffer(final VariantContext vc, final Genotype genotype, final Allele allele, final Integer alleleCode) {
        try {
            alleleBuffer.putInt(alleleCode);
        } catch (final BufferOverflowException e) {
            throw new RuntimeException(
                String.format(
                    "Allele buffer overflow at position: %d code: %d for variant: %s, genotype: %s allele: %s",
                    alleleBuffer.position(),
                    alleleCode,
                    vc.toStringWithoutGenotypes(),
                    genotype == null ? "genotype missing" : genotype.toString(),
                    allele == null ? "no allele present" : allele.toString()),
                e);
        }
    }

    private void updatePhasingBuffer(final VariantContext vc, final Genotype genotype, final byte phaseCode) {
        try {
            phasingBuffer.put(phaseCode);
        } catch (final BufferOverflowException e) {
            throw new RuntimeException(
                String.format(
                    "Phase buffer overflow at position: %d code: %d for variant: %s, genotype: %s",
                    alleleBuffer.position(),
                    phaseCode,
                    vc.toStringWithoutGenotypes(),
                    genotype == null ? "genotype missing" : genotype.toString()),
                e);
        }
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
