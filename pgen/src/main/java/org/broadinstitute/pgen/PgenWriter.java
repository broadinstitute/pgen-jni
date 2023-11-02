/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

package org.broadinstitute.pgen;

import com.github.luben.zstd.ZstdOutputStream;

import htsjdk.io.HtsPath;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
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
     * variant count is not known up front, as long as the corresponding file mode param used is either {@link #PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX}
     * or {@link #PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY} (the write mode must not be {@link #PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK}, which requires an accurate
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
    public static String PVAR_EXTENSION = ".pvar.zst";
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

     /**
     * Enum for representing the subset of plink2 chromosome coding schemes that are supported by this writer.
     * In plink2, the chromosome coding scheme is used when writing various output formats
     * (see https://www.cog-genomics.org/plink/2.0/data#irreg_output).
     * The codes are used here to determine the names of the haploid sex chromosomes in the incoming data set, for purposes
     * of correctly handling haploid calls.
     */
    public enum PgenChromosomeCode {
        /**
         * Autosomes numeric, X/Y single-character, MT two-character, XY/PAR1/PAR2 as usual.
         * Aligns with the b37 reference. This is the default coding used by PLINK 2.
         */
        PLINK_CHROMOSOME_CODE_MT("MT", "X", "Y", "MT"),
        /**
         * PAR1/PAR2 as usual, other chromosomes are 'chr' followed by a numeric code.
         * Aligns with the hg38 reference.
         */
        PLINK_CHROMOSOME_CODE_CHRM("chrM", "chrX", "chrY", "chrM");

        private final String codeString;        // the corresponding chromosome code as recognized by the plink2 command line
        private final String xChromosomeName;   // the chromosome name for the X chromosome for this chromosome code
        private final String yChromosomeName;   // the chromosome name for the Y chromosome for this chromosome code
        private final String mChromosomeName;   // the chromosome name for the mitochondrial chromosome for this chromosome code

        private PgenChromosomeCode(final String codeString, final String xChromosomeName, final String yChromosomeName, final String mChromosomeName) {
             this.codeString = codeString;
             this.xChromosomeName = xChromosomeName;
             this.yChromosomeName = yChromosomeName;
             this.mChromosomeName = mChromosomeName;
        }
        public String value() { return this.codeString; }
        public String getXChromosomeName() { return xChromosomeName; };
        public String getYChromosomeName() { return yChromosomeName; };
        public String getMChromosomeName() { return mChromosomeName; };
    };

    private static final byte PHASED_CODE = (byte) 1;
    private static final byte UNPHASED_CODE = (byte) 0;
    private static final int HAPLOID_PLOIDY = 1;
    private static final int DIPLOID_PLOIDY = 2;

    private final int maxAltAlleles;
    private final boolean lenientPloidyValidation;
    private final List<String> sampleNames;
    private final String xChromosomeName;
    private final String yChromosomeName;
    private final String mChromosomeName;

    private HtsPath pVarFile;
    private HtsPath pSamFile;
    private HtsPath logFile;
    private VariantContextWriter pVarWriter;
    private BufferedWriter logFileWriter;
    private long pgenContextHandle;
    private ByteBuffer alleleBuffer;
    private ByteBuffer phasingBuffer;
    private long expectedVariantCount = 0L;
    private long droppedVariantCount = 0L;
    private long droppedSampleCount = 0L;

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
     * Create a PGEN writer for writing VariantContexts to a Plink2 PGEN file. The writer creates a set of related PGEN files (.pgen and
     * .pvar/.psam files). Depending on the {@code #PgenWriteMode} specified, may also create a .pgen.pgi file.
     * 
     * Supports only diploid sites with fewer than {@link #PLINK2_MAX_ALTERNATE_ALLELES} (this value can be further limited using
     * {@code maxAltAlleles}. Sites that exceed the maximum number of alternate alleles are silently dropped (but will be written to the
     * log file if one is provided). If no log file is provided, any sample that does not have a conforming ploidy (haploid or diploid for
     * sex chromosomes, or stricly diploid for autosomes) will cause an exception to be thrown. Failure to provide a log file does not
     * change the behavior of
     * 
     * @param pgenFileName the name of the PGEN file to be created (must end in .pgen)
     * @param vcfHeader a valid VCF header to use to create the PGEN
     * @param pgenWriteMode the PGEN write mode to use (see {@code PgenWriteMode})
     * @param writeFlags the write flags to use - see {@code PgenWriteFlag}. If phase information is present for the source genotypes, include
     * the {@link PgenWriteFlag#PRESERVE_PHASING} flag. If multi allelic variants are present, include the {@link PgenWriteFlag#MULTI_ALLELIC} flag.
     * @param chromosomeCode the plink2 chromosome coding scheme to use - see {@link PgenChromosomeCode}
     * @param lenientPloidyValidation PGEN requires individual sample to be diploid (except for sex chromsomes, which may be haploid - these are accepted
     * and recoded for pgen as heterozygous/diploid). By default, any ploidy failure will result in an exception to be thrown. Use tru for this value to
     * tolerate ploidy failures (samples will be recoded as missing, and logged if a pg file is provided).
     * @param numberOfVariants the number of variants to be written. if the number is unknown, use the sentinel value {@link PgenWriter#VARIANT_COUNT_UNKNOWN},
     * but doing so precludes the use of the PGEN write mode {@link PgenWriteMode#PGEN_FILE_MODE_BACKWARD_SEEK}
     * @param maxAltAlleles the maximum number of alternate alleles to consider; site with more alterate alleles than this value will be
     * silently dropped
     * @parm logFile any variants or samples that are dropped are logged to this file. may be null.
     **/
    public PgenWriter(
        final HtsPath pgenFileName,
        final VCFHeader vcfHeader,
        final PgenWriteMode pgenWriteMode,
        final EnumSet<PgenWriteFlag> writeFlags,
        final PgenChromosomeCode chromosomeCode,
        final boolean lenientPloidyValidation,
        final long numberOfVariants,
        final int maxAltAlleles,
        final String logFile) {

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

        switch(chromosomeCode) {
            // at the moment, the only difference between the two supported codes is the name of the mitochondrial chromosome, but capture the
            // x and y chromosome names for completeness
            case PLINK_CHROMOSOME_CODE_MT:
            case PLINK_CHROMOSOME_CODE_CHRM:
                this.xChromosomeName = chromosomeCode.getXChromosomeName();
                this.yChromosomeName = chromosomeCode.getYChromosomeName();
                this.mChromosomeName = chromosomeCode.getMChromosomeName();
                break;
            default:
             throw new PgenException(
                String.format("Unrecognized chromosome code name (%s)", chromosomeCode));
        }

        if (logFile != null) {
            this.logFile = new HtsPath(logFile);
            try {
                logFileWriter = Files.newBufferedWriter(this.logFile.toPath());
            } catch (IOException e) {
                throw new RuntimeIOException(String.format("Error opening dropped variants log file %s", logFile), e);
            }
        }
        this.lenientPloidyValidation = lenientPloidyValidation;

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
        
        alleleBuffer = createBuffer(vcfHeader.getNGenotypeSamples() * DIPLOID_PLOIDY * 4); //samples * ploidy * bytes in int32_t (sizeof AlleleCode)
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
        pVarWriter.close();
        pVarWriter = null;

        if (logFileWriter != null) {
            try {
                logFileWriter.close();
            } catch (IOException e) {
                throw new RuntimeIOException(String.format("Error closing dropped variants log file %s", logFile), e);
            }
        }

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
            if (logFileWriter != null) {
                try {
                    logFileWriter.write(String.format("Dropped variant at: %s/%d - too many alleles (%d)\n", vc.getContig(), vc.getStart(), vc.getNAlleles()));
                } catch (IOException e) {
                    throw new RuntimeIOException(String.format("Error writing to dropped variants log file %s", logFile), e);
                }
            }
            return;
        }
        
        alleleBuffer.clear();
        phasingBuffer.clear();
        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);
    
        // Because there may be missing genotyopes, it is significantly simpler to have the primary iteration be
        // through the sample names from the header, rather than through the genotypes.
        // Note: As it stands, this code does not detect or reject the case where there are one or more genotypes
        // in the VC that have a sample name that is not in the header (we could certainly keep track of that and
        // throw, but is it worth the expense ?), or that the order of the genotypes does not match that of the
        // header.
        for (final String sampleName : sampleNames) {
            final Genotype g = vc.getGenotype(sampleName);
            if (g != null) {
                final int ploidy = g.getPloidy();
                if (ploidy == HAPLOID_PLOIDY && (vc.getContig().equals(xChromosomeName) || vc.getContig().equals(yChromosomeName))) {
                    // we have a haploid X or Y, and need to convert it to diploid to satisfy plink
                    final List<Allele> alleles = g.getAlleles();
                    if (alleles.size() != 1) {
                        throw new PgenException(
                            String.format("A genotype with haploid ploidy (%d) does not have one allele (%d) at variant (%s)",
                                ploidy,
                                alleles.size(),
                                vc.toStringWithoutGenotypes()));
                    }
                    final Allele allele = alleles.get(0);
                    final Integer alleleCode = alleleMap.get(allele);
                    if (alleleCode == null) {
                        // do we need this test ? VariantContext doesn't seem to allow such a thing to be created
                        throw new PgenException(
                            String.format("Allele %s not found in allele map for variant %s", allele.toString(), vc.toStringWithoutGenotypes()));
                    }
                    updateAlleleBuffer(vc, g, allele, alleleCode);
                    updateAlleleBuffer(vc, g, allele, alleleCode);
                    updatePhasingBuffer(vc, g, g.isPhased() ? PHASED_CODE : UNPHASED_CODE);
                } else if (ploidy != DIPLOID_PLOIDY) {
                    if (lenientPloidyValidation) {
                        // if lenient, fill in unphased diploid no-call values for any genotype with questionable ploidy
                        updateAlleleBuffer(vc, g, null, PLINK2_NO_CALL_VALUE);
                        updateAlleleBuffer(vc, g, null, PLINK2_NO_CALL_VALUE);
                        updatePhasingBuffer(vc, null, UNPHASED_CODE);
                        if (logFileWriter != null) {
                            try {
                                logFileWriter.write(String.format("Coding non-diploid sample %s as missing at contig/start: %s %d",
                                    g.getSampleName(),
                                    vc.getContig(),
                                    vc.getStart()));
                                droppedSampleCount++;
                            } catch (IOException e) {
                                throw new RuntimeIOException(String.format("Error writing to dropped variants log file %s", logFile), e);
                            }
                        }
                    } else {
                        throw new PgenException(
                            String.format("PGEN only supports diploid calls, but a non-diploid sample (%s) with ploidy (%d) was found at variant (%s)",
                                g.getSampleName(),
                                ploidy,
                                vc.toStringWithoutGenotypes()));
                    }
                } else {
                    for (final Allele allele : g.getAlleles()) {
                        final Integer alleleCode = alleleMap.get(allele);
                        if (alleleCode == null) {
                            // do we need this test ? VariantContext doesn't seem to allow such a thing to be created
                            throw new PgenException(
                                String.format("Allele %s not found in allele map for variant %s", allele.toString(), vc.toStringWithoutGenotypes()));
                        }
                        updateAlleleBuffer(vc, g, allele, alleleCode);
                    }
                    updatePhasingBuffer(vc, g, g.isPhased() ? PHASED_CODE : UNPHASED_CODE);
                }
            } else {
                // fill in unphased diploid no-call values for the missing genotype
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
     * @return the number of variants dropped because they exceeded the max alternate allele count. dropped variants are not written
     * to the .pvar file, but are written to the log file if one was provided.
     */
    public long getDroppedVariantCount() { return droppedVariantCount; }

    /**
     * @return the number of times a sample was dropped because it did not satisfy the ploidy requirements (note that the same sample
     * may be dropped multiple times, once for each site that it fails to satisfy the ploidy requirements for). dropped samples are
     * written to the log file if one was provided.
     */
    public long getDroppedSampleCount() { return droppedSampleCount; }

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
        final OutputStream pVarOutputStream = pVarFile.getOutputStream();
        ZstdOutputStream zstdStream;
        try {
            zstdStream = new ZstdOutputStream(pVarOutputStream);
        } catch (final IOException e) {
            throw new RuntimeIOException(String.format("Error creating the zstd output stream for the .pvar file %s", pVarFile.getRawInputString()), e);
        }
        // technically we're writing a PVAR, not a VCF, but we can do so by composing the VCFWriter with a ZstdOutputStream
        pVarWriter = new VariantContextWriterBuilder()
            .clearOptions()
            .setOptions(EnumSet.of(Options.DO_NOT_WRITE_GENOTYPES, Options.ALLOW_MISSING_FIELDS_IN_HEADER))
            .setOutputStream(zstdStream)
           .build();

        // Ideally there would be a way to record the provenance/origin of a PGEN file right in the file itself, so we can identify
        // files written by this writer, but there isn't. So instead add a "source=..." VCF header line to the .pvar, similar to the
        // "##source=PLINKv2.00" one plink2 adds when it writes a .pvar:
        vcfHeader.addMetaDataLine(new VCFHeaderLine("source", "\"Broad Institute PGEN/PVAR writer\""));
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
