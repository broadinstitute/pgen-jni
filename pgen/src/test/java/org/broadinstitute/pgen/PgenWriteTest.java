/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package org.broadinstitute.pgen;

import htsjdk.io.HtsPath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

import org.broadinstitute.pgen.TestUtils.PgenFileSet;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class PgenWriteTest {

    @SuppressWarnings("unused")
    @Test(expectedExceptions = PgenJniException.class)
    public void testExceptionPropagation() {
        try {
            // force an exception to be thrown from pgen-lib by trying to read from /dev/null
            final PgenWriter unused = new PgenWriter(new HtsPath("/dev/null"), 2, PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES, 6, 3);
        } catch (final PgenJniException e) {
            Assert.assertNotNull(e.getMessage());
            throw e;
        }
    }

    // this test is redundant with the more useful testRoundTripCompareWithPlink2 test below, but is convenient for debugging...
    @Test
    public void testWritePGENSimple() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("test", false);
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                2,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                6,
                3)) {
            for (final VariantContext vc : reader) {
                writer.add(vc);
            }
        }

        // for now, just make sure there are contents
        final long pgenSize = Files.size(pfs.pGenPath());
        Assert.assertNotEquals(pgenSize, 0L);
    }

    @DataProvider(name="roundTripThroughPlink2Tests")
    public Object[][] roundTripThroughPlink2Provider() {
        return new Object[][] {
            // each file mode, without compression
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), 0, false }, // unphased GTs
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), 1, false }, // unphased GTs
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), 2, false }, // unphased GTs

            // this file has phased genotypes, but the concordance test ignores phasing for now since
            // its not preserved by the pgen writer
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), 0, false }, // phased GTs
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), 1, false }, // phased GTs
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), 2, false }  // phased GTs
        };
    }

    @Test(dataProvider = "roundTripThroughPlink2Tests")
    public void testRoundTripCompareWithPlink2(
        final Path originalVCF,
        final int pgenWriteMode,
        final boolean compressPGEN) throws IOException, InterruptedException {
 
        // first, convert the test VCF to pgen twice, once using the PgenWriter and once using plink2
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, compressPGEN);
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF, compressPGEN);

        //TODO: remove this when .pvar/.psam writing are implemented
        // propagate the contents of the plink2-generated .pvar and .psam for now, since generating those isn't implemented
        // yet, and plink2 needs them to do the round-trip back to VCF)
        TestUtils.copyPGENCompanionFiles(plink2FileSet, jniFileSet);

        // now, use plink2 to reconvert both of the pgens back into VCFs
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(jniFileSet, "FromJNI", compressPGEN);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(plink2FileSet, "FromPlink2", compressPGEN);

        if (pgenWriteMode != 1) {
            // while we're at it, run plink2 --pgen-diff on the outputs
            // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupted 
            TestUtils.validatePgen_plink2(jniFileSet);
            TestUtils.validatePgen_plink2(plink2FileSet);
        }

        // also, run a plink2 diff
        TestUtils.pgenDiff_plink2(jniFileSet, jniFileSet);

        // finally, compare the two round tripped vcfs to see if they're equivalent (note that equivalence doesn't mean they're
        // correct,  only that the pgen we generated is concordant with plink2's)
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, vcfFromPGEN_plink2);
    }

    @SuppressWarnings("unused")
    @Test(expectedExceptions = PgenJniException.class)
    public void testInvalidPGENWriteMode() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("invalidPGENwriteMode", false);
        try {
            // don't use try-with-resources here since if we do we'll get a "no variants written" exception instead
            final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                4, // must be one of 0, 1, or 2: 4 is not valid
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES, 
                6,
                3);
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("Invalid pgenWriteMode value"));
            throw e;
         }
    }

    @Test(expectedExceptions = PgenJniException.class)
    public void testClosePGENWithNoWrites() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                1, // use write mode == 1 here so we don't have to clean up the temp file when we abort artificially
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                6, // claim we'll write 6 variants, but don't write them
                3)) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("number of written variants"));
            throw e;
         }
    }

     @Test(expectedExceptions = PgenJniException.class)
    public void testTemporarilyRejectMultiAllelc() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                2,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES, 
                6,
                3)) {
                    final List<Allele> alleles = List.of(
                        Allele.REF_A, Allele.ALT_C, Allele.ALT_G
                    );
                    final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
                    final VariantContext ploidy3VC = vcb.genotypes(
                        List.of(new GenotypeBuilder().name("s1").alleles(alleles).make())
                    ).make();
                    Assert.assertEquals(ploidy3VC.getMaxPloidy(2), 3);
                    pgenWriter.add(ploidy3VC);
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("multi-allelic variants are not yet implemented"));
            throw e;
         }
    }

    // commented out becuase VS Code is too dumb to notice that its disabled
    // // TODO: This test is disabled until we remove the artificial "isBiAllelic" guard that is temporarily in place
    // // to prevent multi-allelics from being processed, since they aren't yet implemented.
    // @Test(expectedExceptions = PgenJniException.class, enabled = false)
    // public void testRejectNonDiploid() throws IOException {
    //     final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
    //     try (final PgenWriter pgenWriter = new PgenWriter(
    //             new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
    //             2,
    //             PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES, 
    //             6,
    //             3)) {
    //                 final List<Allele> alleles = List.of(
    //                     Allele.REF_A, Allele.ALT_C, Allele.ALT_G
    //                 );
    //                 final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
    //                 final VariantContext ploidy3VC = vcb.genotypes(
    //                     List.of(new GenotypeBuilder().name("s1").alleles(alleles).make())
    //                 ).make();
    //                 Assert.assertEquals(ploidy3VC.getMaxPloidy(2), 3);
    //                 pgenWriter.add(ploidy3VC);
    //      } catch (final PgenJniException e) {
    //         Assert.assertTrue(e.getMessage().contains("ploidy = 3"));
    //         throw e;
    //      }
    // }

    @Test(expectedExceptions = PgenJniException.class)
    public void testRequestMaxAltAllelesExceeded() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                2,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES + 1, 
                6,
                3)) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("exceeds the supported pgen max"));
            throw e;
         }
    }

}
