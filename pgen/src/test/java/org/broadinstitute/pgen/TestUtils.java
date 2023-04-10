package org.broadinstitute.pgen;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.testng.Assert;

import htsjdk.io.HtsPath;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

// These utilities require plink2, and assume that it is installed locally and on the path. Since we rely
// on plink2 to verify the files that are written by pgen-jni, it is prefereable to  use a version of
// plink2 that is built from the same source as is used to build pgen-jni, since otherwise in theory
// round-trip errors may surface that are due to the different implementations.

public class TestUtils {
    public static String pGenExtensionUncompressed = ".pgen";
    public static String pGenExtensionCompressed = ".pgen.zst";
    public static String pGenIndexExtension = ".pgen.pgi";    
    public static String pVarExtension = ".pvar";
    public static String pSamExtension = ".psam";
    public static String plinkLogExtension = ".log";
    public static String vcfExtension = FileExtensions.VCF;

    // class/record to hold a set of plink2 companion files as paths (.pgen/.pvar/.psam/.log)
    public record PgenFileSet(Path pGenPath, Path pVarPath, Path pSamPath, Path plinkLogFile) {

        // return the full pathname of the fileset's pgen file, without the file extension
        public String getFileSetPrefix() {
            return getAbsoluteFileNameWithoutExtension(
                pGenPath,
                pGenPath.endsWith(pGenExtensionCompressed) ? pGenExtensionCompressed : pGenExtensionUncompressed);
        }

        // create a trio of temporary pgen (.pgen/.pvar/.psam), plus a log file
        public static PgenFileSet createTempPgenFileSet(final String namePrefix, final boolean useCompressedOutput) throws IOException {
            final String pGenExtension = useCompressedOutput ? pGenExtensionCompressed : pGenExtensionUncompressed;

            final Path pGenPath = createTempFile(namePrefix, pGenExtension).toPath();
            final String pgenNameWithoutExtension  = TestUtils.getLocalFileNameWithoutExtension(pGenPath, pGenExtension);

            // make sure any possible companion files are also marked for deletion
            // we have to force these to be created in order for them to be deleted on exit
            final Path pVarPath = pGenPath.resolveSibling(pgenNameWithoutExtension + pVarExtension).toAbsolutePath();
            pVarPath.toFile().createNewFile();
            pVarPath.toFile().deleteOnExit();

            final Path pSamPath = pGenPath.resolveSibling(pgenNameWithoutExtension + pSamExtension).toAbsolutePath();
            pSamPath.toFile().createNewFile();
            pSamPath.toFile().deleteOnExit();

            final Path pLogPath = pGenPath.resolveSibling(pgenNameWithoutExtension + plinkLogExtension).toAbsolutePath();
            pLogPath.toFile().createNewFile();
            pLogPath.toFile().deleteOnExit();

           // in case we generate an external .pgi index file
           final Path pGenIndexPath = pGenPath.resolveSibling(pgenNameWithoutExtension + pGenIndexExtension).toAbsolutePath();
           pGenIndexPath.toFile().createNewFile();
           pGenIndexPath.toFile().deleteOnExit();

           return new PgenFileSet(pGenPath, pVarPath, pSamPath, pLogPath);
        }
    }
    
    // use pgen-jni to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
    public static PgenFileSet vcfToPgen_jni(
            final Path originalVCF,
            final int pgenFileMode,
            final boolean compressPGEN) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_jni", compressPGEN);
        try(VCFFileReader reader = new VCFFileReader(originalVCF, false)) {
            final int nSamples = reader.getFileHeader().getNGenotypeSamples();
            final List<VariantContext> vcs = new ArrayList<>();
            for (VariantContext vc : reader) {
                vcs.add(vc);
            }
            try (final PgenWriter writer = new PgenWriter(
                                                new HtsPath(pgenFileSet.pGenPath.toAbsolutePath().toString()),
                                                pgenFileMode,
                                                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                                                vcs.size(),
                                                nSamples)) {
                for (final VariantContext vc : vcs) {
                    writer.add(vc);
                }
            }
        }
        return pgenFileSet;
    }

     // use plink2 to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
     public static PgenFileSet vcfToPgen_plink2(final Path originalVCF, final boolean compressPGEN) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_plink2", compressPGEN);
        final String runCommand = String.format("plink2 --vcf %s --make-pgen --out %s", originalVCF.toAbsolutePath(), pgenFileSet.getFileSetPrefix());
 
        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);

        return pgenFileSet;
    }

   // use plink2 to convert a PGEN fileset to a temporary vcf
    public static Path pgenToVCF_plink2(final PgenFileSet pgenFileSet, final String sourceContext, final boolean isCompressed)  throws InterruptedException, IOException {
        final Path plinkGeneratedVCF = createTempFile(pgenFileSet.getFileSetPrefix() + "_" + sourceContext, vcfExtension).toPath();

        // the following plink command will create it's own log file, so make sure it gets marked for deletion as well...
        makeCompanionLogFileTemporary(plinkGeneratedVCF, vcfExtension);

        final String runCommand = String.format(
            "plink2 --pfile %s %s --export vcf --out %s",
            pgenFileSet.getFileSetPrefix(),
            getDecompressCommand(pgenFileSet.pGenPath(), isCompressed),
            TestUtils.getAbsoluteFileNameWithoutExtension(plinkGeneratedVCF, vcfExtension));

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
 
    // Propagate the contents of the .pvar and .psam file from a Plink2-generated fileset to a jni-generated
    // fileset (temporary until the PgenWriter supports writing these companion files).
    //
    // This is a temporary hack needed so we can use plink2 (which requires a complete pgen fileset to convert back
    // to VCF) to validate pgen output without having to generated the .pvar/.psam files ourself. When we have that,
    // this can be removed.
    //
    public static void copyPGENCompanionFiles(final PgenFileSet sourcFileSet, final PgenFileSet targetFileSet) throws IOException {
        // use StandardCopyOption.REPLACE_EXISTING because we've already created the file in order to mark it deleteOnExit
        Files.copy(sourcFileSet.pSamPath, targetFileSet.pSamPath, StandardCopyOption.REPLACE_EXISTING);
        Files.copy(sourcFileSet.pVarPath, targetFileSet.pVarPath, StandardCopyOption.REPLACE_EXISTING);
    }

    // execute an external command
    public static int executeExternalCommand(final String runCommand) throws InterruptedException, IOException {
        final Process process = Runtime.getRuntime().exec(runCommand);
        final int retCode = process.waitFor();
        if (retCode != 0) {
            try (final InputStream errStream = process.getErrorStream();
                 final BufferedReader bufferedReader = new BufferedLineReader(errStream)) {
                while (bufferedReader.ready()) {
                    /// propagate any errors to stderr
                    final String s = bufferedReader.readLine();
                    System.err.println(s);
                }        
            }
        }
        
        return retCode;
    }

    // compare a plink2-generated VCF with a pgen-jni-generated VCF
    public static void verifyRoundTripGenotypeConcordance(final Path plinkVCF, final Path jniVCF) {
        try (final VCFFileReader jniReader = new VCFFileReader(jniVCF, false);
             final CloseableIterator<VariantContext> jniIt = jniReader.iterator()) {
            try (final VCFFileReader plinkReader = new VCFFileReader(plinkVCF, false)) {
                for (final VariantContext plinkVariant : plinkReader) {
                    if (!jniIt.hasNext()) {
                        throw new IllegalStateException("No corresponding jni variant for plink variant: " + plinkVariant);
                    }
                    final VariantContext jniVariant = jniIt.next();
                    assertVariantContextsAreNominallyEqual(plinkVariant, jniVariant);
                }
            }
            Assert.assertFalse(jniIt.hasNext());
        }
    }

    //Asserts that the two provided VariantContext objects have equalsite/position and equal genotypes (other attributes are ignored)
    public static void assertVariantContextsAreNominallyEqual( final VariantContext actual, final VariantContext expected ) {
        Assert.assertNotNull(actual, "VariantContext expected not null");
        Assert.assertEquals(actual.getContig(), expected.getContig(), "chr");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "start");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "end");
        Assert.assertEquals(actual.getID(), expected.getID(), "id");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles for " + expected + " vs " + actual);

        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes(), "hasGenotypes");
        if ( expected.hasGenotypes() ) {
            assertEqualsSet(actual.getSampleNames(), expected.getSampleNames(), "sample names set");
            Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "sample names");
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                assertGenotypesAreEqual(actual.getGenotype(sample), expected.getGenotype(sample));
            }
        }
    }

    public static final <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new HashSet<T>(actual);
        final Set<T> expectedSet = new HashSet<T>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    // Asserts that the two provided Genotype objects are concordant (not identical, only the same name, alleles and type)
    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype sample names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        //phasing differences are caught here...
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");
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

    // given a Path, return the absolute path of the file, without the trailing extension
    public static String getAbsoluteFileNameWithoutExtension(final Path targetPath, final String extension) {
        final String targetAbsolutePath = targetPath.toAbsolutePath().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf(extension));
    }

     // given a Path, return the prefix of only the terminal part of the filename (only the last segment of the
     // filename, without any extension)
     public static String getLocalFileNameWithoutExtension(final Path targetPath, final String extension) {
        final String targetAbsolutePath = targetPath.getFileName().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf(extension));
    }
 
    //TODO: the decompress arg doesn't seem to be recognized by the current plink2 code
    private static String getDecompressCommand(final Path pgenPath, final boolean isCompressed) {
        return isCompressed ?
            "--zst-decompress " +  pgenPath.toAbsolutePath().toString() + " ":
            " ";
    }

    // generate a potential plink2-generated logfile name make sure its marked for deletion
    private static Path makeCompanionLogFileTemporary(final Path plink2File, final String plink2FileExtension) {
        final File conversionLogFile = new File(plink2File.resolveSibling(
            TestUtils.getLocalFileNameWithoutExtension(plink2File, plink2FileExtension) + plinkLogExtension).toFile().toString());
        conversionLogFile.deleteOnExit();
        return conversionLogFile.toPath();
    }

}