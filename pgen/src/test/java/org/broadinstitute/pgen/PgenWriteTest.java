/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

package org.broadinstitute.pgen;

import htsjdk.io.HtsPath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

import org.broadinstitute.pgen.PgenWriter.PgenWriteFlag;
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
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

public class PgenWriteTest {

    @SuppressWarnings("unused")
    @Test(expectedExceptions = PgenException.class)
    public void testExceptionPropagation() throws IOException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("pgenReadOnly");
        pgenFileSet.pGenPath().toFile().setReadOnly();
        try {
            // force an exception to be thrown from pgen-lib by trying to write to a file that is read only
            final PgenWriter unused = new PgenWriter(
                new HtsPath(pgenFileSet.getFileSetPrefix()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES);
        } catch (final PgenException e) {
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
                    EnumSet.noneOf(PgenWriteFlag.class),
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
            // Cases for a test that creates a PGEN from a VCF, and then validates the generated PGEN by:
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

            // format
            // file name, write mode, use true variant count, extra plink args, write flags
            
            // small, all bi-allelic, unphased (6 variants/3 samples), test once for each write mode, and for write mode != PGEN_FILE_MODE_BACKWARD_SEEK,
            // also test with variant count provided up front or not (PGEN_FILE_MODE_BACKWARD_SEEK always requires an acccurate variant count)
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, true, "--output-chr M", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr M", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr M", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, true, "--output-chr M", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, false, "--output-chr M", EnumSet.noneOf(PgenWriteFlag.class) },

            // slightly larger, all bi-allelic, phased (600 variants/2504 samples); the genotype concordance validation ignores
            // phasing for now since its not preserved by the pgen writer)
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, true, "--output-chr M", EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr M", EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr M", EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, true, "--output-chr M", EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, false, "--output-chr M", EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
 
            // larger still, unphased, includes ~6000 multiallelic sites (~117,932 variants/10 samples)
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },

            {  Paths.get("testdata/hg38_trio.pik3ca.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            {  Paths.get("testdata/hg38_trio.pik3ca.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            // same as hg38_trio.pik3ca.vcf above, but with the genotypes for the (multi-allelic) variant at 179135392 modified so the last site allele
            // is not referenced by any genotype, and the (multi-allelic) variant at site 179170076 modified so the middle (index 1) allele is not
            // referenced by any genotype; this triggers the issue described https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ
            { Paths.get("testdata/hg38_trio.pik3ca.unreferenced.allele.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, true, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/hg38_trio.pik3ca.unreferenced.allele.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
 
            // TODO: add multi-allelic AND phased tests

            // temporary/local test cases
            // { Paths.get("testdata/external/0000000009-NHGRI_AnVIL_3K.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            // { Paths.get("testdata/external/0000000004-NHGRI_AnVIL_3K.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) },
            // { Paths.get("testdata/external/0000000002-NHGRI_AnVIL_3K.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, false, "--output-chr chr26", EnumSet.noneOf(PgenWriteFlag.class) }
        };
    }

    @Test(dataProvider = "roundTripCompareWithPlink2Provider")
    public void testRoundTripCompareWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final boolean useTrueVariantCount,
        final String extraPlinkArgs,
        final EnumSet<PgenWriteFlag> writeFlags) throws IOException, InterruptedException {

        // the pgenlib api can't handle PGEN_FILE_MODE_BACKWARD_SEEK unless a true variant count is provided, so catch
        // attempts to test this combination
        Assert.assertTrue(useTrueVariantCount == true || pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK);
 
        // first, convert the test VCF to pgen twice; once using the PgenWriter and once using plink2 directly
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, useTrueVariantCount, writeFlags);
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
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, vcfFromPGEN_plink2, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING));

        // finally, for extra measure, compare the pgen-jni round-tripped VCF with the ORIGINAL VCF
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, originalVCF, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING));
        System.out.println("Done with: " + originalVCF);
    }

    // ensure the PgenWriter constructor rejects attempts to use VARIANT_COUNT_UNKNOWN with PGEN_FILE_MODE_BACKWARD_SEEK file mode,
    // since its forbidden by plink2
    @Test(expectedExceptions = PgenException.class)
    public void testSeekWriteModeRejectsUnknownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testSeekWriteModeRejectsUnknownVariantCount");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));

        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES))
        {
            // do nothing...
        } catch (final PgenException e) {
            Assert.assertTrue(e.getMessage().contains("requires a known variant count"));
            throw e;
        }
    }
 
    @Test
    public void testNoWritesWithKnownVariantCount() throws IOException {
        // this test is basically to ensure that the pgen-lib C++ code correctly handles closing in the case where
        // the variant count is known, but 0 variants are written (the plink2 code seems to handle it ok if too few
        // variants are written, but not if 0 variants are written, so this tests that pgen-lib specifically handles
        // that case
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testNoWritesWithKnownVariantCount");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                EnumSet.noneOf(PgenWriteFlag.class),
                6, // claim we'll write 6 variants, but don't write them
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES))
        {
            // do nothing...
            Assert.assertEquals(pgenWriter.getDroppedVariantCount(), 0L);
            Assert.assertEquals(pgenWriter.getWrittenVariantCount(), 0L);
        }
    }

    @Test
    public void testTooFewWritesWithKnownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testTooFewWritesWithKnownVariantCount");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                EnumSet.noneOf(PgenWriteFlag.class),
                vcfMetaData.nVariants(),
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES))
        {
            // write 1 variant, then bail
            pgenWriter.add(reader.iterator().next());
            Assert.assertEquals(pgenWriter.getDroppedVariantCount(), 0L);
            Assert.assertEquals(pgenWriter.getWrittenVariantCount(), 1L);
        }
    }

    // Create a writer with VARIANT_COUNT_UNKNOWN (so it will not throw on close if you haven't written enough variants),
    // and then never write any variants. This should succeed with file mode PGEN_FILE_MODE_WRITE_SEPARATE_INDEX (it will
    // fail with filemode PGEN_FILE_MODE_BACKWARD_SEEK, since in that case you're required to provide a variant count up
    // front, and then write that many variants).
    @Test(expectedExceptions=PgenEmptyPgenException.class)
    public void testNoWritesWithUnknownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testNoWritesWithUnknownVariantCount");
        try (final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    TestUtils.createSingleSampleVCFHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    EnumSet.noneOf(PgenWriteFlag.class),
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES);
            final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false))
        {
            // do nothing
        }
     }

   @DataProvider(name="plink2CapabilitiesProvider")
    public Object[][] plink2CapabilitiesProvider() {
        return new Object[][] {
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26", true },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chr26", true },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr M", false },
        };
    }

    // disabled, since this is not really a pgen-jni test, just a convenient way to see how plink2 handles conversions
    @Test(enabled = false, dataProvider = "plink2CapabilitiesProvider")
    public void testPlink2Capabilities(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final String additionalPlinkArgs,
        final boolean ignorePhasing) throws IOException, InterruptedException {
        // convert the test VCF to pgen, then back to vcf using plink2
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(plink2FileSet, "FromPlink2", additionalPlinkArgs);

        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
            // while we're at it, run plink2 --pgen-diff on the output
            // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupted 
            TestUtils.validatePgen_plink2(plink2FileSet);
        }

        // compare the round tripped vcf with the original to see if they're concordant
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_plink2, originalVCF, ignorePhasing);
    }

    // reject non-diploid variants
    @Test(expectedExceptions = PgenException.class)
    public void testRejectNonDiploid() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("rejectNonDiploid");
        final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
        final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
        final VariantContext ploidy3VC = vcb.genotypes(
            List.of(new GenotypeBuilder().name("s1").alleles(alleles).make())
        ).make();
        Assert.assertEquals(ploidy3VC.getMaxPloidy(2), 3);

        // For reasons unknown, putting the PgenWriter constructor inside the try-with-resources results
        // in an IllegalArgumentException, due to an illegal attempt to self-suppress an exception
        // (somewhere code is trying to add this exception to it's own suppressed exception list), but only
        // when runnning inside VS Code. This does not happen when gradle runs the test, but keep it
        // outside.
        PgenWriter pgenWriter =  null;
        try {
            pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES);
        } catch (final PgenException e) {
            // we expect a PgenException in this test, but not from the constructor, so make sure we don't
            // inadvertently succeed because of some nefarious failure in the constructor
            Assert.fail("should never reach here");
        }
        try {
            pgenWriter.add(ploidy3VC);
        } catch (final PgenException ploidyException) {
            Assert.assertTrue(ploidyException.getMessage().contains("ploidy = 3"));
            try {
                pgenWriter.close();
            } catch (final PgenEmptyPgenException emptyPgenException) {
                // this is expected because the ploidy exception is thrown on the very first variantwe try to
                // add, so no variants ever get written
            }
            throw ploidyException;
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
                EnumSet.noneOf(PgenWriteFlag.class),
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

    // verify that we reject attempts to set a max alternate allele threshold that exceeds the plink2 maximum
    @Test(expectedExceptions = PgenException.class)
    public void testRejectRequestedMaxAltAllelesExceedsPlinkMax() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectRequestedMaxAltAllelesExceedsPlinkMax");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                6,
                // add one to plink's max to be certain we exceed the limit
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES + 1)) {
         } catch (final PgenException e) {
            Assert.assertTrue(e.getMessage().contains("exceeds the supported PGEN max"));
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
                    vcfMetaData.vcfHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    EnumSet.noneOf(PgenWriteFlag.class),
                    // use VARIANT_COUNT_UNKNOWN, since otherwise we'll get an exception because we didn't write
                    // enough variants when the try-with-resources calls close on the writer
                    PgenWriter.VARIANT_COUNT_UNKNOWN, 
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES)) {
                // first, add one good variant
                VariantContext vc = reader.iterator().next();
                writer.add(vc);

                // now conjure and write up a bad variant that will cause buffer overflow due to having too many genotypes
                final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                final List<Genotype> gts = new ArrayList<>();
                vc.getGenotypes().stream().forEach(gt -> gts.add(gt));
                vc.getGenotypes().stream().forEach(gt -> gts.add(gt));

                writer.add(vcb.genotypes(gts).make()); // throws a runtime exeption due to buffer overflow exception
        } catch (final RuntimeException e) {
            // The underlying writer code wraps the BufferOverflowException in a RuntimeException that is decorated
            // with additional context information. Make sure this RuntimeException is actually caused by a
            // BufferOverflowException, and then rethrow.
            Assert.assertTrue(e.getCause() instanceof BufferOverflowException);
            Assert.assertTrue(e.getMessage().contains("Buffer overflow at position:"));
            throw e;
        }
    }

}
