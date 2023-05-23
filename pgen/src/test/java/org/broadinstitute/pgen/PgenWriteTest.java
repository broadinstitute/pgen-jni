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
import java.nio.BufferOverflowException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class PgenWriteTest {

    @SuppressWarnings("unused")
    @Test(expectedExceptions = PgenJniException.class)
    public void testExceptionPropagation() throws IOException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("pgenReadOnly");
        pgenFileSet.pGenPath().toFile().setReadOnly();
        try {
            // force an exception to be thrown from pgen-lib by trying to write to a file that is read only
            final PgenWriter unused = new PgenWriter(
                new HtsPath(pgenFileSet.getFileSetPrefix()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                6,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES);
        } catch (final PgenJniException e) {
            Assert.assertNotNull(e.getMessage().contains("kPglRetOpenFail"));
            throw e;
        }
    }

    // this test is redundant with the more useful testRoundTripCompareWithPlink2 test below, but is convenient for debugging...
    @Test
    public void testWriteSimpleBiallelic() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("writeBiallelic");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    vcfMetaData.vcfHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    vcfMetaData.nVariants(),
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES
                    )) {
                reader.forEach(vc -> writer.add(vc));
            }

        // for now, just make sure there are contents
        final long pgenSize = Files.size(pfs.pGenPath());
        Assert.assertNotEquals(pgenSize, 0L);
    }

    @DataProvider(name="roundTripCompareWithPlink2Provider")
    public Object[][] roundTripCompareWithPlink2Provider() {
        return new Object[][] {
            // This test creates a PGEN from a VCF, and then validates the generated PGEN by:
            //
            //  1) Using plink2 to directly generate a PGEN from the same input VCF.
            //  2) Using plink2 to directly validate both the plink2-generated and pgen-jni-generated PGENs.
            //  3) Using plink2 to diff the pgen-jni-generated PGEN with the plink2-generated PGEN, verifying
            //     that no differences are reported.
            //  4) Using plink2 to regenerate VCFs from the two generated PGENs, and then:
            //      a) running a concordance test to validate that the VCF generated from the pgen-jni-generated PGEN is equivalent
            //         to the VCF generated from the plink2-generated PGEN
            //      b) running a concordance test to validate that the pgen-jni-generated VCF is equivalent to the original test VCF
            //
            // Since plink2 doesn't respect the chromosome names contained in the PGEN's companion PVAR when generating VCFs,
            // and instead generates names using one of several predefined schemes identified by codes that can be provided
            // on the command line via the "--output-chr" argument, each of these test cases has to include an appropriate
            // "--output-chr" argument in order to make the subsequent VCF comparison to the original succeed. See
            // https://www.cog-genomics.org/plink/2.0/data#irreg_output.
            
            // small, all bi-allelic, unphased (6 variants/3 samples), test once for each write mode, all without compression
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, "--output-chr M" },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr M" },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, "--output-chr M" },

            // slightly larger, all bi-allelic, phased (600 variants/2504 samples); the genotype concordance validation ignores
            // phasing for now since its not preserved by the pgen writer)
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, "--output-chr M" },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr M" },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, "--output-chr M" },
 
            // larger still, includes ~6000 multiallelic sites (~117,932 variants/10 samples)
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26" },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26" },

            {  Paths.get("testdata/hg38_trio.pik3ca.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26" },
            // same as hg38_trio.pik3ca.vcf above, but with the genotypes for the variant at 179135392 modified so the last site allele is not
            // referenced by any genotype; this triggers the issue described https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ
            // currently this fails because we haven't picked up the new plink code with the fix (described here
            // https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ); this should be fixed the next time we update to new code
            // and a newer version of plink2, and then this test can be enabled
            { Paths.get("testdata/hg38_trio.pik3ca.unreferenced.allele.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26" }
        };
    }

    @Test(dataProvider = "roundTripCompareWithPlink2Provider")
    public void testRoundTripCompareWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final String extraPlinkArgs) throws IOException, InterruptedException {
 
        // first, convert the test VCF to pgen twice; once using the PgenWriter and once using plink2 directly
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode);
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF);

        // now use plink2 to validate both of the generated pgen file sets
        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
            // while we're at it, run plink2 --pgen-diff on the outputs
            // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupt 
            TestUtils.validatePgen_plink2(jniFileSet);
            TestUtils.validatePgen_plink2(plink2FileSet);
        }

        // run a plink2 diff on the two filesets and make sure there are no differences reported
        TestUtils.pgenDiff_plink2(jniFileSet, plink2FileSet);

        // now, use plink2 to reconvert both of the generated pgens back into VCFs, and then compare the two round-tripped VCFs to
        // each other see if they're concordant (note that concordance doesn't mean they're correct, only that the VCF created from
        // the pgen fileset that was created by the pgen writer is concordant with the VCF created from the pgen fileset created by
        // plink2)
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(jniFileSet, "FromJNI", extraPlinkArgs);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(plink2FileSet, "FromPlink2", extraPlinkArgs);
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, vcfFromPGEN_plink2);

        // finally, for extra measure, compare the pgen-jni round-tripped VCF with the ORIGINAL VCF
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, originalVCF);
    }

   @Test(expectedExceptions = PgenJniException.class)
    public void testClosePGENWithNoWrites() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("noWritesPgenTest");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                 // use PGEN_FILE_MODE_WRITE_AND_COPY mode here so we don't have to clean up the temp file when we abort artificially
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                6, // claim we'll write 6 variants, but don't write them
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES
            )) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("number of variants written"));
            throw e;
         }
    }

    // reject non-diploid variants
    @Test(expectedExceptions = PgenJniException.class)
    public void testRejectNonDiploid() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("rejectNonDiploid");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                6,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES)) {
            final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
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

    // verify that we really do drop variants that exceed the maximum alterate allele threshold
    @Test
    public void testRejectTooManyAltAlleles() throws IOException {
        final int ARTIFICALLY_LOW_MAX_ALLELE_THRESHOLD = 3;
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("tooManyAltAlelesTest");
 
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                vcfMetaData.nVariants() + 1, // add one to account for the synthetic variant that we add that gets dropped
                ARTIFICALLY_LOW_MAX_ALLELE_THRESHOLD)) {

            try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false)) {
                reader.forEach(vc -> pgenWriter.add(vc));
            }

            Assert.assertEquals(pgenWriter.getWrittenVariantCount(), vcfMetaData.nVariants());
            Assert.assertEquals(pgenWriter.getDroppedVariantCount(), 0);

            // now write a variant that should get dropped due to exceeding the artifically low max allele threshold we set above
            final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G, Allele.ALT_T);
            final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
            final VariantContext tooManyAlleles = vcb.genotypes(List.of(new GenotypeBuilder().name("s1").alleles(alleles).make())).make();
            pgenWriter.add(tooManyAlleles);

            Assert.assertEquals(pgenWriter.getWrittenVariantCount(), vcfMetaData.nVariants());
            Assert.assertEquals(pgenWriter.getDroppedVariantCount(), 1);    
         }
    }

    // verify that we reject attempts to set a max alternate allele threshold that exceeds plink2 maximum
    @Test(expectedExceptions = PgenJniException.class)
    public void testRejectRequestedMaxAltAllelesExceedsPlinkMax() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("maxAltAlelesTest");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                6,
                // add one to plink's max to be certain we exceed the limit
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES + 1)) {
         } catch (final PgenJniException e) {
            Assert.assertTrue(e.getMessage().contains("exceeds the supported pgen max"));
            throw e;
         }
    }

    // force our internal allele buffer (the one containing the allele_codes that are passed to plink2) to be overflowed
    @Test(expectedExceptions = RuntimeException.class)
    public void testAlleleBufferOverflow() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("writeBufferOverflow");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    // use our test header, which has only a single sample (this will cause an artifically small allele_code buffer
                    // to be allocated, which will then be overflowed when using the variants from the 3-sample test file)
                    TestUtils.createSingleSampleVCFHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    vcfMetaData.nVariants(),
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES)) {
            reader.forEach(vc -> writer.add(vc));
        } catch (final RuntimeException e) {
            // The underlying writer code wraps the BufferOverflowException in a RuntimeException that is
            // decorated with additional context information. Make sure this RuntimeException is actually caused by a
            // BufferOverflowException, and then rethrow.
            Assert.assertTrue(e.getCause() instanceof BufferOverflowException);
            Assert.assertTrue(e.getCause().getMessage().contains("Buffer overflow at position:"));
            throw e;
        }
    }

}
