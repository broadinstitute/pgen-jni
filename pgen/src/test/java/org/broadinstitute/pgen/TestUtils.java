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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TestUtils {

    // record to hold a trio of plink files as paths (.pgen/.pvar/.psam)
    public record PgenFileSet(Path pGenPath, Path pVarPath, Path pSamPath) {

        // return the full pathname of the fileset's pgen file, without the file extension
        public String getFileSetPrefix() {
            return getAbsoluteFileNameWithoutExtension(pGenPath);
        }

        // create a trio of temporary pgen (.pgen/.pvar/.psam), plus a log file
        public static PgenFileSet createTempPgenFileSet(final String namePrefix) throws IOException {
            final Path pGenPath = createTempFile(namePrefix, ".pgen").toPath();
            
            final String pgenNameWithoutExtension  = TestUtils.getLocalFileNameWithoutExtension(pGenPath);

            final Path pVarPath = pGenPath.resolveSibling(pgenNameWithoutExtension + ".pvar").toAbsolutePath();
            pVarPath.toFile().createNewFile();
            pVarPath.toFile().deleteOnExit();

            // for these, we have to force them to be created in order for them to be deleted on exit
            final Path pSamPath = pGenPath.resolveSibling(pgenNameWithoutExtension + ".psam").toAbsolutePath();
            pSamPath.toFile().createNewFile();
            pSamPath.toFile().deleteOnExit();

            final Path pLogPath = pGenPath.resolveSibling(pgenNameWithoutExtension + ".log").toAbsolutePath();
            pLogPath.toFile().createNewFile();
            pLogPath.toFile().deleteOnExit();

            return new PgenFileSet(pGenPath, pVarPath, pSamPath);
        }
    }
    
    // use pgen-jni to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
    public static PgenFileSet vcfToPgen_jni(final Path originalVCF, final int pgenFileMode) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_jni");
        try(VCFFileReader reader = new VCFFileReader(originalVCF, false)) {
            final int nSamples = reader.getFileHeader().getNGenotypeSamples();
            final List<VariantContext> vcs = new ArrayList<>();
            for (VariantContext vc : reader) {
                vcs.add(vc);
            }
            try (final PgenWriter writer = new PgenWriter(
                                                new HtsPath(pgenFileSet.pGenPath.toAbsolutePath().toString()),
                                                pgenFileMode,
                                                vcs.size(),
                                                nSamples)) {
                for (VariantContext vc : vcs) {
                    writer.add(vc);
                }
            }
        }

        return pgenFileSet;
    }

     // use plink2 to convert a VCF to PGEN and return the resulting (temporary) files as a PgenFileSet
     public static PgenFileSet vcfToPgen_plink2(final Path originalVCF) throws IOException, InterruptedException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("vcfToPgen_plink2");
        final String runCommand = String.format("plink2 --vcf %s --make-pgen --out %s", originalVCF.toAbsolutePath(), pgenFileSet.getFileSetPrefix());
 
        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);

        return pgenFileSet;
    }

   // use plink2 to convert a PGEN fileset to a temporary vcf
    public static Path pgenToVCF_plink2(final PgenFileSet pgenFileSet, final String sourceContext)  throws InterruptedException, IOException {
        final Path plinkGeneratedVCF = createTempFile(pgenFileSet.getFileSetPrefix() + "_" + sourceContext, ".vcf").toPath();

        // the following plink command will create it's own log file, so make sure that gets marked for deletion as well...
        final File conversionLogFile = new File(plinkGeneratedVCF.resolveSibling(
            TestUtils.getLocalFileNameWithoutExtension(plinkGeneratedVCF) + ".log").toFile().toString());
        conversionLogFile.deleteOnExit();

        final String runCommand = String.format(
            "plink2 --pfile %s --export vcf --out %s",
            pgenFileSet.getFileSetPrefix(),
            TestUtils.getAbsoluteFileNameWithoutExtension(plinkGeneratedVCF));

        final int cmdResult = executeExternalCommand(runCommand);
        Assert.assertEquals(cmdResult, 0);

        return plinkGeneratedVCF;
    }
      
    // Propagate the contents of the .pvar and .psam file from a Plink2-generated fileset to a jni-generated
    // fileset (temporary until the PgenWriter supports writing these companion files).
    //
    // This is a temporary hack needed so we can use plink2 (which requires a complete pgen fileset to convert back
    // to VCF) to validate pgen output without having to generated the .pvar/.psam files ourself. When we have that,
    // this can be removed.
    //
    public static void propagatePlink2FileContents(final PgenFileSet plink2FileSet, final PgenFileSet jniFileSet) throws IOException {
        // use StandardCopyOption.REPLACE_EXISTING because we've already created the file in order to mark it deleteOnExit
        Files.copy(plink2FileSet.pSamPath, jniFileSet.pSamPath, StandardCopyOption.REPLACE_EXISTING);
        Files.copy(plink2FileSet.pVarPath, jniFileSet.pVarPath, StandardCopyOption.REPLACE_EXISTING);
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

    // Creates a temp file that will be deleted on exit after tests are complete.
    public static File createTempFile(final String name, final String extension) {
        try {
            File file = File.createTempFile(name, extension);
            file.deleteOnExit();
            return file;
        } catch (IOException ex) {
            throw new RuntimeException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    // given a Path, return the prefix of the full (absolute) filename (the relative filename without any extension)
    public static String getAbsoluteFileNameWithoutExtension(final Path targetPath) {
        final String targetAbsolutePath = targetPath.toAbsolutePath().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf("."));
    }

    // given a Path, return the prefix of the terminal filename (only the last segment of the filename without any extension)
    public static String getLocalFileNameWithoutExtension(final Path targetPath) {
        final String targetAbsolutePath = targetPath.getFileName().toString();
        return targetAbsolutePath.substring(0, targetAbsolutePath.lastIndexOf("."));
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

    //TODO: does this corectly handle phase ?
    // Asserts that the two provided Genotype objects are concordant (not identical, only the same name, alleles and type)
    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");
    }

}