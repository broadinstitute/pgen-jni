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

import org.broadinstitute.pgen.PgenWriter.PgenWriteMode;
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
        final File readOnlyFile = TestUtils.createTempFile("pgenReadOnly", ".pgen");
        readOnlyFile.setReadOnly();
        try {
            // force an exception to be thrown from pgen-lib by trying to write to a file that is read only
            final PgenWriter unused = new PgenWriter(
                new HtsPath(readOnlyFile.getAbsolutePath()),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                6,
                3);
        } catch (final PgenJniException e) {
            Assert.assertNotNull(e.getMessage().contains("kPglRetOpenFail"));
            throw e;
        }
    }

    // this test is redundant with the more useful testRoundTripCompareWithPlink2 test below, but is convenient for debugging...
    @Test
    public void testWritePGENBiallelic() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("test", false);
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                    vcfMetaData.nVariants(),
                    vcfMetaData.nSamples())) {
                reader.forEach(vc -> writer.add(vc));
            }

        // for now, just make sure there are contents
        final long pgenSize = Files.size(pfs.pGenPath());
        Assert.assertNotEquals(pgenSize, 0L);
    }

    @DataProvider(name="roundTripThroughPlink2Tests")
    public Object[][] roundTripThroughPlink2Provider() {
        return new Object[][] {
            // small, bi-allelic, unphased - once for each file mode, without compression
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, false },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, false },

            // slightly larger, bi-allelic, phased (the concordance validation ignores phasing for
            // now since its not preserved by the pgen writer)
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, false },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, false },
 
            // has a few multi-allelic variants - not used until multi-allelic varints are implemented
            //{ Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false },
            //{ Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false },
       };
    }

    @Test(dataProvider = "roundTripThroughPlink2Tests")
    public void testRoundTripCompareWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final boolean compressPGEN) throws IOException, InterruptedException {
 
        // first, convert the test VCF to pgen twice, once using the PgenWriter and once using plink2
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, compressPGEN);
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF, compressPGEN);

        //TODO: remove this when .pvar/.psam writing are implemented
        // propagate the contents of the plink2-generated .pvar and .psam for now, since generating those isn't implemented
        // yet, and plink2 needs them to do the round-trip back to VCF)
        TestUtils.copyPGENCompanionFiles(plink2FileSet, jniFileSet);

        // now, use plink2 to reconvert both of the pgens back into VCFs
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(jniFileSet, "FromJNI", null, compressPGEN);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(plink2FileSet, "FromPlink2", null, compressPGEN);

        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
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

    @DataProvider(name="plink2CapabilitiesProvider")
    public Object[][] plink2CapabilitiesProvider() {
        return new Object[][] {
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY },
        };
    }

    // disabled, since this is not really a test, just a convenient way to see how plink2 handles conversions
    @Test(dataProvider = "plink2CapabilitiesProvider", enabled = false)
    public void testPlink2Capabilities(final Path originalVCF, final PgenWriteMode pgenWriteMode) throws IOException, InterruptedException {
        // convert the test VCF to pgen, then back to vcf using plink2
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF, false);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(
            plink2FileSet, "FromPlink2",
             "--output-chr chr26",      // tell plink how to name contigs in the output
              false);

        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
            // while we're at it, run plink2 --pgen-diff on the output
            // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupted 
            TestUtils.validatePgen_plink2(plink2FileSet);
        }

        // compare the round tripped vcf with the original to see if they're concordant
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_plink2, originalVCF);
    }

   @Test(expectedExceptions = PgenJniException.class)
    public void testClosePGENWithNoWrites() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                 // use write and copy mode here so we don't have to clean up the temp file when we abort artificially
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES,
                6, // claim we'll write 6 variants, but don't write them
                3)) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("number of written variants"));
            throw e;
         }
    }

     @Test(expectedExceptions = PgenJniException.class)
    public void testTemporarilyRejectMultiAllelic() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
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

    // TODO: This test is disabled until we remove the artificial "isBiAllelic" guard that is temporarily in place
    // to prevent multi-allelics from being processed, since they aren't yet implemented.
    @Test(expectedExceptions = PgenJniException.class, enabled = false)
    public void testRejectNonDiploid() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
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
            Assert.assertTrue(e.getMessage().contains("ploidy = 3"));
            throw e;
         }
    }

    @Test(expectedExceptions = PgenJniException.class)
    public void testRequestMaxAltAllelesExceeded() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest", false);
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                PgenWriter.MAX_PLINK2_ALTERNATE_ALLELES + 1, 
                6,
                3)) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("exceeds the supported pgen max"));
            throw e;
         }
    }

}
