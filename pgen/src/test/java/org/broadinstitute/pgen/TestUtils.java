/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

package org.broadinstitute.pgen;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.pgen.PgenWriter.PgenChromosomeCode;
import org.broadinstitute.pgen.PgenWriter.PgenWriteMode;
import org.testng.Assert;

import htsjdk.io.HtsPath;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

// These utilities require plink2, and assume that it is installed locally and on the path. Since we rely
// on plink2 to verify the files that are written by pgen-jni, it is prefereable to  use a version of
// plink2 that is built from the same source as is used to build pgen-jni, since otherwise in theory
// round-trip errors may surface that are due to the different implementations.

public class TestUtils {

    public static String plinkLogExtension = ".log";
    public static final String SINGLE_SAMPLE_HEADER_SAMPLE_NAME = "atLeastOneSampleRequired";

    // class/record to hold a set of temporary plink2 companion files as paths (.pgen/.pvar/.psam/.log)
    public record PgenFileSet(Path pGenPath, Path pVarPath, Path pSamPath, Path plinkLogFile) {

        // create a trio of temporaary pgen files (.pgen/.pvar/.psam), and mark them, and any other possible companion files, for deletion
        public static PgenFileSet createTempPgenFileSet(final String namePrefix) throws IOException {
            final String pGenExtension = PgenWriter.PGEN_EXTENSION;

            final Path pGenPath = createTempFile(namePrefix, pGenExtension).toPath();
            final String pgenNameWithoutExtension  = getLocalFileNameWithoutExtension(pGenPath, pGenExtension);

            // we have to force these to be created in order for them to be deleted on exit
            final Path pVarPath = pGenPath.resolveSibling(pgenNameWithoutExtension + PgenWriter.PVAR_EXTENSION).toAbsolutePath();
            pVarPath.toFile().createNewFile();
            pVarPath.toFile().deleteOnExit();

            final Path pSamPath = pGenPath.resolveSibling(pgenNameWithoutExtension + PgenWriter.PSAM_EXTENSION).toAbsolutePath();
            pSamPath.toFile().createNewFile();
            pSamPath.toFile().deleteOnExit();

            // make sure any other possible companion files are also marked for deletion, even if they aren't used

            // the .log file
            final Path pLogPath = pGenPath.resolveSibling(pgenNameWithoutExtension + plinkLogExtension).toAbsolutePath();
            pLogPath.toFile().createNewFile();
            pLogPath.toFile().deleteOnExit();

           // the .pgi index file
           final Path pGenIndexPath = pGenPath.resolveSibling(pgenNameWithoutExtension + PgenWriter.PGEN_INDEX_EXTENSION).toAbsolutePath();
           pGenIndexPath.toFile().createNewFile();
           pGenIndexPath.toFile().deleteOnExit();

           // the .pgen.tmp file (in case the file mode used in the test results in one)
           final Path pGenTempPath= pGenPath.resolveSibling(pgenNameWithoutExtension + PgenWriter.PGEN_EXTENSION + ".tmp").toAbsolutePath();
           pGenTempPath.toFile().createNewFile();
           pGenTempPath.toFile().deleteOnExit();

           return new PgenFileSet(pGenPath, pVarPath, pSamPath, pLogPath);
        }

        // return the full pathname of the fileset's pgen file, without the file extension
        public String getFileSetPrefix() { return PgenWriter.getAbsoluteFileNameWithoutExtension(pGenPath, PgenWriter.PGEN_EXTENSION); }

    }

    // record to hold metadata for a VCF
    public record VcfMetaData(VCFHeader vcfHeader, long nVariants){};

    @SuppressWarnings("unused")
    public static VcfMetaData getVcfMetaData(final Path vcfPath) {
        long nVariants = 0L;
        int nSamples = 0;
        VCFHeader vcfHeader;
        try (final VCFFileReader vcfReader = new VCFFileReader(vcfPath, false);
             final CloseableIterator<VariantContext> vcIt = vcfReader.iterator()) {
             nVariants = vcIt.stream().count();
            vcfHeader = vcfReader.getFileHeader();
        }
        return new VcfMetaData(vcfHeader, nVariants);
    }

    // create a VCFHeader that has at least one sample, since pgen will assert (in a debug build) if the number of samples is 0
    @SuppressWarnings("unchecked")
    public static VCFHeader createSingleSampleVCFHeader() {
        final VCFHeader vcfHeader = new VCFHeader(
            Collections.EMPTY_SET,
            new HashSet<String>() {{
                add(SINGLE_SAMPLE_HEADER_SAMPLE_NAME);
            }});
        Assert.assertEquals(vcfHeader.getNGenotypeSamples(), 1, "the test ehaders must have at least one sample");
        return vcfHeader;
    }
    
    // use pgen-jni to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
    public static PgenFileSet vcfToPgen_jni(
            final Path originalVCF,
            final PgenWriteMode pgenWriteMode,
            final PgenChromosomeCode chromosomeCode,
            final boolean useTrueVariantCount,
            final EnumSet<PgenWriter.PgenWriteFlag> writeFlags) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_jni");
        final VcfMetaData vcfMetaData = getVcfMetaData(originalVCF);
        try(final VCFFileReader reader = new VCFFileReader(originalVCF, false);
            final PgenWriter writer = new PgenWriter(
                    new HtsPath(pgenFileSet.pGenPath.toAbsolutePath().toString()),
                    vcfMetaData.vcfHeader,
                    pgenWriteMode,
                    writeFlags,
                    chromosomeCode,
                    false,
                    useTrueVariantCount == true ? vcfMetaData.nVariants : PgenWriter.VARIANT_COUNT_UNKNOWN,
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                    null)) {
            reader.forEach(vc -> writer.add(vc));
            // display the variant counts
            System.out.println(String.format(
                    "%d variants written\n%d variants dropped",
                    writer.getWrittenVariantCount(),
                    writer.getDroppedVariantCount()));
        }
        return pgenFileSet;
    }

     // use plink2 to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
     public static PgenFileSet vcfToPgen_plink2(final Path originalVCF, final PgenChromosomeCode chromosomeCode) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_plink2");
        final String runCommand = String.format(
            "plink2 --vcf %s --output-chr %s --make-pgen --out %s",
            originalVCF.toAbsolutePath(),
            chromosomeCode.value(),
            pgenFileSet.getFileSetPrefix());
 
        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);

        return pgenFileSet;
    }

   // Use plink2 to convert a PGEN fileset to a temporary vcf
    public static Path pgenToVCF_plink2(
        final PgenFileSet pgenFileSet,
        final String sourceContext,
        final String additionalArgs) throws InterruptedException, IOException {
        final Path plinkGeneratedVCF = createTempFile(pgenFileSet.getFileSetPrefix() + "_" + sourceContext, FileExtensions.VCF).toPath();

        // the following plink command will create it's own log file, so make sure it gets marked for deletion as well...
        makeCompanionLogFileTemporary(plinkGeneratedVCF, FileExtensions.VCF);
        
        final String runCommand = String.format(
            "plink2 --pfile %s %s --export vcf --out %s",
            pgenFileSet.getFileSetPrefix(),
            additionalArgs == null ? "" : additionalArgs,
            PgenWriter.getAbsoluteFileNameWithoutExtension(plinkGeneratedVCF, FileExtensions.VCF));

        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);

        return plinkGeneratedVCF;
    }

    // Use plink2 to validate a pgen:
	//     plink2 --pfile pgenFile --validate
    public static void validatePgen_plink2(final PgenFileSet pgenFileSet) throws IOException, InterruptedException {
        final String runCommand = String.format(
            "plink2 --pfile %s --validate",
            pgenFileSet.getFileSetPrefix());

        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);
    }

    // Use plink2 to generate a Pgen diff:
	//     plink2 --pfile pgenFile1 --pgen-diff pgenFile2
    public static void pgenDiff_plink2(final PgenFileSet sourceFileSet, final PgenFileSet targetFileSet) throws IOException, InterruptedException {
        final String runCommand = String.format(
            "plink2 --pfile %s --pgen-diff %s",
            sourceFileSet.getFileSetPrefix(),
            targetFileSet.getFileSetPrefix());

        // Unlike other commands, which seem to generate a log file adajcent to the command inputs, with a name derived
        // from the inputs, --pgen-diff generates a plink2.log file and a plink2.pdiff file in the current directory, so
        // make sure they're marked for deletion
        final File plink2_pDiffFile = new File("plink2.pdiff");
        plink2_pDiffFile.deleteOnExit();
        final File plink2_diffLogFile = new File("plink2.log");
        plink2_diffLogFile.deleteOnExit();

        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0); 
        
        // verify that there is only 1 line (the header line) in the plink2.pdiff file; otherwise there are
        // actual differences that need to be investigated
        final List<String> diffOutput = IOUtil.slurpLines(plink2_pDiffFile);
        Assert.assertEquals(diffOutput.size(), 1);
        Assert.assertTrue(diffOutput.get(0).equals("#ID\tIID\tGT1\tGT2"));
    }
 
    // execute an external command
    public static int executeExternalCommand(final String runCommand) throws InterruptedException, IOException {
        final Process process = Runtime.getRuntime().exec(runCommand);
        final int retCode = process.waitFor();
        if (retCode != 0) {
            // propagate any errors to stderr
            try (final InputStream errStream = process.getErrorStream();
                 final BufferedReader errBufferedReader = new BufferedLineReader(errStream)) {
                while (errBufferedReader.ready()) {
                    System.err.println(errBufferedReader.readLine());
                }        
            }
        }
        
        return retCode;
    }

    // compare two VCFs for genotype concordance - note that this does not compare the VCFHeaders
    public static void verifyRoundTripGenotypeConcordance(final Path actualVCF, final Path expectedVCF, final boolean ignorePhasing, final boolean ignoreGenotypeSwaps) {        
        long phaseMatches = 0L;
        long concordantPhaseFailures = 0L;
        try (final VCFFileReader expReader = new VCFFileReader(expectedVCF, false);
             final CloseableIterator<VariantContext> expIt = expReader.iterator()) {
            try (final VCFFileReader actReader = new VCFFileReader(actualVCF, false)) {
                final VCFHeader expHeader = expReader.getFileHeader();
                final VCFHeader actHeader = actReader.getFileHeader();
                final HashMap<String, Integer> actOffsetMap = actHeader.getSampleNameToOffset();
                final HashMap<String, Integer> expOffsetMap = expHeader.getSampleNameToOffset();

                // ensure that our samples are in the same order in the headers
                for (final Map.Entry<String, Integer> entry : actOffsetMap.entrySet()) {
                    Assert.assertEquals(entry.getValue(), expOffsetMap.get(entry.getKey()));
                }

                for (final VariantContext actVariant : actReader) {
                    if (!expIt.hasNext()) {
                        throw new IllegalStateException("No corresponding jni variant for plink variant: " + actVariant);
                    }
                    final VariantContext expVariant = expIt.next();
                    final long phaseFailures = assertVariantContextsAreNominallyEqual(actVariant, expVariant, ignorePhasing, ignoreGenotypeSwaps);
                    concordantPhaseFailures += phaseFailures;
                    phaseMatches += (actHeader.getNGenotypeSamples() - phaseFailures); 
               }
            }
            Assert.assertFalse(expIt.hasNext());
        }
        if (!ignorePhasing) {
            // when phasing is turned on, plink2 roundtrips any unphased genotypes that are either homozygous or no-call as phased, so report
            // those as harmless (if there are any other phase mismatches, then the concordance check, and thus the whole test, will fail)
            System.out.println(String.format("Phase matches: %d Homozygous/no-call phase changes introduced by plink2: %d", phaseMatches, concordantPhaseFailures));
        }
    }

    //Asserts that the two provided VariantContext objects have equal site/position and concordant genotypes (other attributes are ignored)
    public static long assertVariantContextsAreNominallyEqual( final VariantContext actual, final VariantContext expected, final boolean ignorePhasing, final boolean ignoreGenotypeSwaps ) {
        Assert.assertNotNull(actual, "VariantContext expected not null");
        Assert.assertEquals(actual.getContig(), expected.getContig(), "chr");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "start");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "end");
        Assert.assertEquals(actual.getID(), expected.getID(), "id");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles for " + expected + " vs " + actual);

        long concordantPhaseFailures = 0L;
        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes(), "hasGenotypes");
        if ( expected.hasGenotypes() ) {
            assertEqualsSet(actual.getSampleNames(), expected.getSampleNames(), "sample names set");
            Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "sample names");
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                concordantPhaseFailures += assertGenotypesAreConcordant(actual.getGenotype(sample), expected.getGenotype(sample), ignorePhasing, ignoreGenotypeSwaps);
            }
        }
        return concordantPhaseFailures;
    }

    public static final <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new HashSet<T>(actual);
        final Set<T> expectedSet = new HashSet<T>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    // Asserts that the two provided Genotype objects are concordant (not necessarily identical, only that the alleles and type match,
    // and that phasing matches for called variants (because plink2 changes homozygous variants to be phased, even if they're not, so
    // we count these just to keep track, but accept them as concordant).
    public static long assertGenotypesAreConcordant(final Genotype actual, final Genotype expected, final boolean ignorePhasing, final boolean ignoreGenotypeSwaps) {
        // if you're going to ignore genotype swaps, you better also ignore phasing...
        Assert.assertTrue(ignoreGenotypeSwaps == false || ignorePhasing == true, "bad test argument combination");

        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype sample names");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        long concordantPhaseFailures = 0L;
        if (ignoreGenotypeSwaps) {
            // just check that the same alleles are present, even if swapped
            assertEqualsSet(new HashSet<>(actual.getAlleles()), new HashSet<>(expected.getAlleles()), "Genotype alleles");
        } else {
            Assert.assertTrue(expected.sameGenotype(actual, ignorePhasing), "Genotype alleles");
            if (!ignorePhasing) {
                final String actGenString = actual.getGenotypeString();
                final String expGenString = expected.getGenotypeString();
                if (!actGenString.equals(expGenString)) {
                    // plink2 roundtrips unphased genotypes that are either homozygous, or  no-calls, as phased, whether the input
                    // says they were phased or not, so detect and track those, but accept them as "concordant"
                    Assert.assertEquals(actual.getAllele(0), expected.getAllele(0));
                    Assert.assertEquals(actual.getAllele(1), expected.getAllele(1));
                    // now make sure they're homozygous (even if they're no-calls), since we only want to accept phase changes as
                    // "concordant" if they're homozygous. note that we test the allele directly since htsjdk doesn't return
                    // isHom()==true for "./." or ".|.", since they're neither homref nor homvar.
                    Assert.assertEquals(actual.getAllele(0), actual.getAllele(1));
                    concordantPhaseFailures++;
                } else {
                    Assert.assertEquals(actGenString, expGenString);
                }
            }
        }

        return concordantPhaseFailures;
    }

     // Creates a temp file that will be deleted on exit after tests are complete.
     public static File createTempFile(final String name, final String extension) {
        try {
            final File file = File.createTempFile(name, extension);
            file.deleteOnExit();
            return file;
        } catch (final IOException ex) {
            throw new RuntimeException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

     // given a Path, return the prefix of only the terminal part of the filename (only the last segment of the
     // filename, without any extension)
     public static String getLocalFileNameWithoutExtension(final Path targetPath, final String extension) {
        final String targetAbsolutePath = targetPath.getFileName().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf(extension));
    }
 
    // generate a potential plink2-generated logfile name make sure its marked for deletion
    private static Path makeCompanionLogFileTemporary(final Path plink2File, final String plink2FileExtension) {
        final File conversionLogFile = new File(plink2File.resolveSibling(
            getLocalFileNameWithoutExtension(plink2File, plink2FileExtension) + plinkLogExtension).toFile().toString());
        conversionLogFile.deleteOnExit();
        return conversionLogFile.toPath();
    }

}