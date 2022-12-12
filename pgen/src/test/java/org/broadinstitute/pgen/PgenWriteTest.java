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

import org.broadinstitute.pgen.PgenWriter.PgenChromosomeCode;
import org.broadinstitute.pgen.PgenWriter.PgenWriteFlag;
import org.broadinstitute.pgen.PgenWriter.PgenWriteMode;
import org.broadinstitute.pgen.TestUtils.PgenFileSet;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

public class PgenWriteTest {

    @SuppressWarnings("unused")
    @Test(expectedExceptions = PgenException.class)
    public void testExceptionPropagation() throws IOException {
        final PgenFileSet pgenFileSet = PgenFileSet.createTempPgenFileSet("testExceptionPropagation");
        pgenFileSet.pGenPath().toFile().setReadOnly();
        try {
            // force an exception to be thrown from pgen-lib by trying to write to a file that is read only
            final PgenWriter unused = new PgenWriter(
                new HtsPath(pgenFileSet.getFileSetPrefix()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,  // doesn't matter
                false,
                PgenWriter.VARIANT_COUNT_UNKNOWN,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                null);
        } catch (final PgenException e) {
            Assert.assertNotNull(e.getMessage().contains("kPglRetOpenFail"));
            throw e;
        }
    }

    // this test is redundant with the more useful testRoundTripCompareWithPlink2 test below, but is convenient for debugging...
    @Test
    public void testWriteSimpleBiallelic() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testWriteSimpleBiallelic");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    vcfMetaData.vcfHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    EnumSet.noneOf(PgenWriteFlag.class),
                    PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,
                    false,
                    vcfMetaData.nVariants(),
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                    null
                    )) {
                reader.forEach(vc -> writer.add(vc));
            }

        // for now, just make sure there are contents
        final long pgenSize = Files.size(pfs.pGenPath());
        Assert.assertNotEquals(pgenSize, 0L);
    }

    @DataProvider(name="roundTripAutosomesWithPlink2Provider")
    public Object[][] roundTripAutosomesWithPlink2Provider() {
        return new Object[][] {
            // Cases for autosome tests that create a PGEN from a VCF, and then validate the generated PGEN by:
            //
            //  1) Using plink2 to directly generate a PGEN from the same input VCF.
            //  2) Using plink2 to directly validate both the plink2-generated and pgen-jni-generated PGENs.
            //  3) Using plink2 to diff the pgen-jni-generated PGEN with the plink2-generated PGEN, verifying
            //     that no differences are reported.
            //  4) Using plink2 to regenerate VCFs from the two generated PGENs, and then:
            //      a) running a concordance test to validate that the VCF generated from the pgen-jni-generated PGEN is concordant
            //          with the VCF generated from the plink2-generated PGEN
            //      b) running a concordance test to validate that the pgen-jni-generated VCF is concordant with the original test VCF
            //
            // Since plink2 doesn't respect the chromosome names contained in the PGEN's companion PVAR when generating VCFs,
            // and instead generates names using one of several predefined schemes identified by codes that can be provided
            // on the command line via the "--output-chr" argument, each of these test cases has to include an appropriate value
            // for the "--output-chr" argument in order to make the subsequent VCF comparison to the original succeed. See
            // https://www.cog-genomics.org/plink/2.0/data#irreg_output.

            // Note that this provide is for autosomes only - sex chromosomes require special handling and are in separate providers.

            // provider structure: file name, write mode, chromosome code, use true variant count, write flags
            
            // small, all bi-allelic, unphased (6 variants/3 samples), test once for each write mode, and for write mode != PGEN_FILE_MODE_BACKWARD_SEEK,
            // also test with variant count provided up front or not (PGEN_FILE_MODE_BACKWARD_SEEK always requires an acccurate variant count)
            // uses the b37 reference, so use the --output-chr code "MT"; technically this is the plink2 default, but secify it anyway just for the record
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, false, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/CEUtrioTest.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, false, EnumSet.noneOf(PgenWriteFlag.class) },

            // slightly larger, all bi-allelic, phased (600 variants/2504 samples); the genotype concordance validation ignores
            // phasing for now since its not preserved by the pgen writer)
            // uses the b37 reference, so use the --output-chr code "MT"; technically this is the plink2 default, but secify it anyway just for the record
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true, EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true,EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, false, EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, true, EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, false, EnumSet.of(PgenWriteFlag.PRESERVE_PHASING) },
 
            // larger still, unphased, includes ~6000 multiallelic sites (~117,932 variants/10 samples)
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) },

            {  Paths.get("testdata/hg38_trio.pik3ca.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, true, EnumSet.noneOf(PgenWriteFlag.class) },
            {  Paths.get("testdata/hg38_trio.pik3ca.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) },
            // same as hg38_trio.pik3ca.vcf above, but with the genotypes for the (multi-allelic) variant at 179135392 modified so the last site allele
            // is not referenced by any genotype, and the (multi-allelic) variant at site 179170076 modified so the middle (index 1) allele is not
            // referenced by any genotype; this triggers the issue described https://groups.google.com/g/plink2-users/c/Sn5qVCyDlDw/m/GOWScY6tAQAJ
            { Paths.get("testdata/hg38_trio.pik3ca.unreferenced.allele.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, true, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/hg38_trio.pik3ca.unreferenced.allele.vcf").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) },
 
            // // multiallelic, partially phased (phasing synthesized by randomly mutating the genotypes in 0000000000-my_demo_filters.vcf.gz)
            { Paths.get("testdata/0000000000-my_demo_filters.partiallyphased.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.of(PgenWriteFlag.MULTI_ALLELIC, PgenWriteFlag.PRESERVE_PHASING) },
        };
    }

    @Test(dataProvider = "roundTripAutosomesWithPlink2Provider")
    public void testRoundTripAutosomesWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final PgenChromosomeCode chromosomeCode,
        final boolean useTrueVariantCount,
        final EnumSet<PgenWriteFlag> writeFlags) throws IOException, InterruptedException {
        final String extraPlinkArgs = "--output-chr " + chromosomeCode.value();

        // the pgenlib api can't handle PGEN_FILE_MODE_BACKWARD_SEEK unless a true variant count is provided, so catch
        // attempts to test this combination
        Assert.assertTrue(useTrueVariantCount == true || pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK);
 
        // first, convert the test VCF to pgen twice; once using the PgenWriter and once using plink2 directly
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, chromosomeCode, useTrueVariantCount, writeFlags);
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF, chromosomeCode);

        // while we're at it, run plink2 --validate on the outputs
        // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupt 
        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
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
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, vcfFromPGEN_plink2, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING), false);

        // finally, for extra measure, compare the pgen-jni round-tripped VCF with the ORIGINAL VCF
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, originalVCF, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING), false);
    }

    @DataProvider(name="roundTripXChromosomeWithPlink2Provider")
    public Object[][] roundTripXChromosomeWithPlink2Provider() {
        return new Object[][] {
            // X chromosome tests; we need to treat validation of variants on the X chromosome a little differently from the tests for autosomal
            // or Y chromosomes because plink2 will not create a vcf from a pgen/pvar if they contain variants on chrX unless the sex for each
            // sample is provided as part of the .psam.
            //
            // Additionally, our X chromosome test data contains genotypes that are written in the order "2/1", which plink2 will rewrite as "1/2",
            // so we need to use a special concordance test that tolerates this. And, since genotypes written that way require that we ignore
            // genotype swaps, we also need to ignore phasing, which in turn requires us to not attempt to use the multi-allelic flag (since for
            // some reason plink only accepts that flag if the phasing is preserved).
            { Paths.get("testdata/shard_109_chrX.vcf.gz").toAbsolutePath(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) },
            { Paths.get("testdata/shard_110_chrX.vcf.gz").toAbsolutePath(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.noneOf(PgenWriteFlag.class) }
        };
    }

    @Test(dataProvider = "roundTripXChromosomeWithPlink2Provider")
    public void testRoundTripXChromosomeWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final PgenChromosomeCode chromosomeCode,
        final boolean useTrueVariantCount,
        final EnumSet<PgenWriteFlag> writeFlags) throws IOException, InterruptedException {
        final String extraPlinkArgs = "--output-chr " + chromosomeCode.value();

        // the pgenlib api can't handle PGEN_FILE_MODE_BACKWARD_SEEK unless a true variant count is provided, so catch
        // attempts to test this combination
        Assert.assertTrue(useTrueVariantCount == true || pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK);
 
        // first, convert the test VCF to pgen using the PgenWriter
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, chromosomeCode, useTrueVariantCount, writeFlags);

        // now, use plink2 to reconvert the generated pgen back into a VCF, and then compare the round-tripped VCF to
        // the original to see if they're concordant (note that concordance doesn't mean it's correct, only that the VCF created from
        // the PGEN fileset that was created by the PGEN writer is concordant with the orignal VCF)
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(jniFileSet, "FromJNI", extraPlinkArgs);

        // Do a special concordance test that is specific to the X chromosome test, because the test data contains some genotypes that are
        // written as "2/1", which plink2 rewrites as "1/2", so use a concordance test that tolerates this.
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, originalVCF, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING), true);
    }

    @DataProvider(name="roundTripYChromosomeWithPlink2Provider")
    public Object[][] roundTripYChromosomeWithPlink2() {
        return new Object[][] {
            // Y chromosome tests; we need to treat these a little differently than the tests that use autosomal or X chromosomes because
            // plink2 will not diff two PGENs if they contain calls on chrY unless you provide the sex for each sample in a .psam.
            //
            // Also, note that although we can do a concordance test between the two intermediate vcfs that come from using plink to
            // roundtrip the two PGENs back to vcf, we cannot do the final concordance test against the original VCF for these test cases,
            // because if the original has chrY encoded as diploid, plink2 will generate vcfs with chrY as haploid. Rather than trying to
            // customize the test infrastructure to accommodate this, we just don't do the final concordance test, and just rely on the
            // concordance between the two (pgenWriter and plink2) conversion paths.

            // Test case with diploid chrY calls
            { Paths.get("testdata/shard_111_chrY.vcf.gz").toAbsolutePath(),
                 PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.of(PgenWriteFlag.MULTI_ALLELIC, PgenWriteFlag.PRESERVE_PHASING) },

            // Test case with haploid chrY calls to test the PGEN writer's "haploid->homozygous diploid" conversion.
            // The test file was created by converting shard_111_chrY.vcf.gz to PGEN with plink2, and then back to .vcf also with plink2
            // (which has the effect of converting all of the "diploid" chrY calls to haploid, even though they are encoded in the PGEN
            // as diploid - plink2 knows they're aligned to a Y chromosome so it apparently converts them back to haploid in the vcf).
            // Also, since plink writes a version 4.3 vcf, I manually changed the version number in this file back to v4.2 so that we could
            // write it back out as a .pvar using htsjdk.
            { Paths.get("testdata/shard_111_chrY_plinked.vcf").toAbsolutePath(),
                 PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, PgenChromosomeCode.PLINK_CHROMOSOME_CODE_CHRM, false, EnumSet.of(PgenWriteFlag.MULTI_ALLELIC, PgenWriteFlag.PRESERVE_PHASING) },
        };
    }

    @Test(dataProvider = "roundTripYChromosomeWithPlink2Provider")
    public void testRoundTripYChromosomeWithPlink2(
        final Path originalVCF,
        final PgenWriteMode pgenWriteMode,
        final PgenChromosomeCode chromosomeCode,
        final boolean useTrueVariantCount,
        final EnumSet<PgenWriteFlag> writeFlags) throws IOException, InterruptedException {
        final String extraPlinkArgs = "--output-chr " + chromosomeCode.value();

        // the pgenlib api can't handle PGEN_FILE_MODE_BACKWARD_SEEK unless a true variant count is provided, so catch
        // attempts to test this combination
        Assert.assertTrue(useTrueVariantCount == true || pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK);
 
        // first, convert the test VCF to pgen twice; once using the PgenWriter and once using plink2 directly
        final TestUtils.PgenFileSet jniFileSet = TestUtils.vcfToPgen_jni(originalVCF, pgenWriteMode, chromosomeCode, useTrueVariantCount, writeFlags);
        final TestUtils.PgenFileSet plink2FileSet = TestUtils.vcfToPgen_plink2(originalVCF, chromosomeCode);

        // while we're at it, run plink2 --validate on the outputs
        // but only if write mode != kPgenWriteSeparateIndex, since that causes plink2 to say the pgen file is corrupt 
        if (pgenWriteMode != PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX) {
            TestUtils.validatePgen_plink2(jniFileSet);
            TestUtils.validatePgen_plink2(plink2FileSet);
        }

        // ...note that for the Y chromosome tests, we skip the plink2 diff on the two filesets, since plink won't diff them unless we
        // provide sex info for each sample

        // now, use plink2 to reconvert both of the generated PGENs back into VCFs, and then compare the two round-tripped VCFs to
        // each other see if they're concordant (note that concordance doesn't necessrily mean they're correct, only that the one
        // created by plink from the pgen file we wrote is concordant with the one created from the pgen that plink wrote).
        //
        // Also note that for the Y chromosome, we don't do a concordance test with the ORIGINAL file (only between the two generated
        // files) because plink2 will write Y chromosome genotypes to a VCF as haploid, even if they were encoded as diploid in the original
        // file (and despite the fact that plink requires us to write them as homozygous/diplod in the pgen), so a concordance test against
        // the original file would fail if it has diploid genotypes, which some of these test inputs have.
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(jniFileSet, "FromJNI", extraPlinkArgs);
        final Path vcfFromPGEN_plink2 = TestUtils.pgenToVCF_plink2(plink2FileSet, "FromPlink2", extraPlinkArgs);
        TestUtils.verifyRoundTripGenotypeConcordance(vcfFromPGEN_jni, vcfFromPGEN_plink2, !writeFlags.contains(PgenWriteFlag.PRESERVE_PHASING), false);
    }

    // ensure the PgenWriter constructor rejects attempts to use VARIANT_COUNT_UNKNOWN with PGEN_FILE_MODE_BACKWARD_SEEK file mode,
    // since its forbidden by plink2
    @Test(expectedExceptions = PgenException.class)
    public void testRejectSeekWriteModeWithUnknownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectSeeWriteModeWithUnknownVariantCount");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));

        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_BACKWARD_SEEK,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,  // doesn't matter
                false,
                PgenWriter.VARIANT_COUNT_UNKNOWN,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                null))
        {
            // do nothing...
        } catch (final PgenException e) {
            Assert.assertTrue(e.getMessage().contains("requires a known variant count"));
            throw e;
        }
    }

    @Test
    public void testAcceptNoWritesWithKnownVariantCount() throws IOException {
        // this test is basically to ensure that the pgen-lib C++ code correctly handles closing in the case where
        // the variant count is known, but 0 variants are written (the plink2 code seems to handle it ok if too few
        // variants are written, but not if 0 variants are written, so this tests that pgen-lib specifically handles
        // that case)
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testAcceptNoWritesWithKnownVariantCount");
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,  // doesn't matter
                false,
                6, // claim we'll write 6 variants, but don't write them
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                null))
        {
            // do nothing...
            Assert.assertEquals(pgenWriter.getDroppedVariantCount(), 0L);
            Assert.assertEquals(pgenWriter.getWrittenVariantCount(), 0L);
        }
    }

    @Test
    public void testAcceptTooFewWritesWithKnownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testAcceptTooFewWritesWithKnownVariantCount");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,
                false,
                vcfMetaData.nVariants(),
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                null))
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
    public void testRejectNoWritesWithUnknownVariantCount() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectNoWritesWithUnknownVariantCount");
        try (final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    TestUtils.createSingleSampleVCFHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    EnumSet.noneOf(PgenWriteFlag.class),
                    PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, // doesn't matter
                    false,
                    PgenWriter.VARIANT_COUNT_UNKNOWN,
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                    null);
            final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false))
        {
            // do nothing
        }
     }

   @DataProvider(name="plink2CapabilitiesProvider")
    public Object[][] plink2CapabilitiesProvider() {
        return new Object[][] {
            { Paths.get("testdata/0000000000-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chrM", true },
            { Paths.get("testdata/0000000001-my_demo_filters.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr chrM", true },
            { Paths.get("testdata/1kg_phase3_chr21_start.vcf.gz").toAbsolutePath(), PgenWriteMode.PGEN_FILE_MODE_WRITE_AND_COPY, "--output-chr M", false },
        };
    }

    // reject non-diploid samples if lenientPloidyValidation==false
    @Test(expectedExceptions = PgenException.class)
    public void testNonDiploidStrict() throws IOException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectNonDiploid");
        final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
        final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
        final VariantContext ploidy3VC = vcb.genotypes(
            List.of(new GenotypeBuilder().name(TestUtils.SINGLE_SAMPLE_HEADER_SAMPLE_NAME).alleles(alleles).make())
        ).make();
        Assert.assertEquals(ploidy3VC.getMaxPloidy(2), 3);

        // For reasons unknown, putting the PgenWriter constructor inside the try-with-resources results
        // in an IllegalArgumentException, but only when runnning inside VS Code, due to an illegal attempt
        // to self-suppress an exception (somewhere code is trying to add this exception to it's own suppressed
        // exception list). This does not happen when gradle runs the test, but keep it outside so the test
        // works in VSCode.
        PgenWriter pgenWriter =  null;
        try {
            pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT, // doesn't matter
                false,
                PgenWriter.VARIANT_COUNT_UNKNOWN,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                null);
        } catch (final PgenException e) {
            // we expect a PgenException in this test, but not from the constructor, so make sure we don't
            // inadvertently succeed because of some nefarious failure in the constructor
            Assert.fail("should never reach here");
        }
        try {
            pgenWriter.add(ploidy3VC);
        } catch (final PgenException ploidyException) {
            Assert.assertTrue(ploidyException.getMessage().contains("ploidy (3)"));
            try {
                pgenWriter.close();
            } catch (final PgenEmptyPgenException emptyPgenException) {
                // this is expected because the ploidy exception is thrown on the very first variant we try to
                // add, so no variants ever get written
            }
            throw ploidyException;
        }
    }

    // verify that if lenientPloidyValiation==true, we mark non-diploid samples as missing just for the one site; and verify 
    // that we count how many of these we drop, and that the event is written to a log file if one is provided
    @Test
    public void testNonDiploidLenient() throws IOException, InterruptedException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectNonDiploid");
        final List<Allele> alleles = List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G);
        final VariantContextBuilder vcb = new VariantContextBuilder("test", "chr1", 1, 1, alleles);
        final VariantContext ploidy3VC = vcb.genotypes(
            List.of(new GenotypeBuilder().name(TestUtils.SINGLE_SAMPLE_HEADER_SAMPLE_NAME).alleles(alleles).make())
        ).make();
        Assert.assertEquals(ploidy3VC.getMaxPloidy(2), 3);

        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                TestUtils.createSingleSampleVCFHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,
                true, // LENIENT required for this test
                PgenWriter.VARIANT_COUNT_UNKNOWN,
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                pfs.pgenLogPath().toAbsolutePath().toString())) {
            Assert.assertEquals(pgenWriter.getDroppedSampleCount(), 0);
            pgenWriter.add(ploidy3VC);
            Assert.assertEquals(pgenWriter.getDroppedSampleCount(), 1);
        }

        // now verify that we logged the dropped variant to the log file
        final String logString = TestUtils.readTextFile(pfs.pgenLogPath());
        Assert.assertTrue(logString.contains("Coding non-diploid sample atLeastOneSampleRequired as missing at contig/start: chr1 1"));

        // now convert the pgen back to VCF and verify that the sample is now marked as missing
        final String extraPlinkArgs = "--output-chr " + PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT.value();
        final Path vcfFromPGEN_jni = TestUtils.pgenToVCF_plink2(pfs, "FromJNI", extraPlinkArgs);

        try (final VCFFileReader vcfReader = new VCFFileReader(vcfFromPGEN_jni, false)) {
            for (final VariantContext v : vcfReader) {
                Assert.assertTrue(v.getGenotype(TestUtils.SINGLE_SAMPLE_HEADER_SAMPLE_NAME).isNoCall());
            }
        }
    }

    // verify that we drop variants that exceed the maximum alterate allele threshold, that we count how many
    // of these we drop, and that we write them to a log file if one is provided
    @Test
    public void testDropSitesWithTooManyAltAlleles() throws IOException {
        final int ARTIFICALLY_LOW_MAX_ALLELE_THRESHOLD = 3;
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testRejectTooManyAltAlleles");
 
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final PgenWriter pgenWriter = new PgenWriter(
                new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                vcfMetaData.vcfHeader(),
                PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                EnumSet.noneOf(PgenWriteFlag.class),
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,    // doesn't matter
                false,
                vcfMetaData.nVariants() + 1, // add one to account for the synthetic variant that we add that gets dropped
                ARTIFICALLY_LOW_MAX_ALLELE_THRESHOLD,
                pfs.pgenLogPath().toAbsolutePath().toString())) {

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
        
        // now verify that we logged the dropped variant to the log file
        final String logString = TestUtils.readTextFile(pfs.pgenLogPath());
        Assert.assertTrue(logString.contains("Dropped variant at: chr1/1 - too many alleles (4)"));
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
                PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,    // doesn't matter
                false,
                6,
                // add one to plink's max to be certain we exceed the limit
                PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES + 1,
                null)) {
         } catch (final PgenException e) {
            Assert.assertTrue(e.getMessage().contains("exceeds the supported PGEN max"));
            throw e;
         }
    }

    @DataProvider(name="missingGenotypes")
    public Object[][] missingGenotypesProvider() {
        return new Object[][] {
            // genotype index for genotypes that should be missing
            { List.of(0) },
            { List.of(1) },
            { List.of(2) },
            { List.of(0, 1) },
            { List.of(1, 2) },
            { List.of(0, 2) },
            { List.of(0, 1, 2) },
        };
    }

    @Test(dataProvider = "missingGenotypes")
    public void testMissingGenotypes(final List<Integer> requestedMissingGenotypeIndices) throws IOException, InterruptedException {
        final PgenFileSet pfs = PgenFileSet.createTempPgenFileSet("testMissingGenotypes");
        final TestUtils.VcfMetaData vcfMetaData = TestUtils.getVcfMetaData(Paths.get("testdata/CEUtrioTest.vcf"));
        try (final VCFFileReader reader = new VCFFileReader(new File("testdata/CEUtrioTest.vcf"), false);
             final PgenWriter writer = new PgenWriter(
                    new HtsPath(pfs.pGenPath().toAbsolutePath().toString()),
                    vcfMetaData.vcfHeader(),
                    PgenWriteMode.PGEN_FILE_MODE_WRITE_SEPARATE_INDEX,
                    EnumSet.noneOf(PgenWriteFlag.class),
                    PgenChromosomeCode.PLINK_CHROMOSOME_CODE_MT,
                    false,
                    PgenWriter.VARIANT_COUNT_UNKNOWN, 
                    PgenWriter.PLINK2_MAX_ALTERNATE_ALLELES,
                    null)) {
                        
            for (final VariantContext vc : reader) {
                // write the variant with missing genotypes to the PGEN
                final VariantContextBuilder vcb = new VariantContextBuilder(vc);
                final List<Genotype> oldGenotypes = new ArrayList<>(vc.getGenotypes());
                final List<Genotype> remainingGenotypes = new ArrayList<>();
                final Map<String, Integer> sampleNameToOffset = vcfMetaData.vcfHeader().getSampleNameToOffset();
                for (final Genotype g : oldGenotypes) {
                    if (!requestedMissingGenotypeIndices.contains(sampleNameToOffset.get(g.getSampleName()))) {
                        remainingGenotypes.add(oldGenotypes.get(sampleNameToOffset.get(g.getSampleName())));
                    }
                }
                writer.add(vcb.genotypes(remainingGenotypes).make());
            }
        }

        // now, round-trip the pgen back to a VCF, and verify that all of the missing genotypes are no-call
        final Path plinkGeneratedVCF = TestUtils.pgenToVCF_plink2(pfs, "missing_genotypes", "--output-chr M");
        try (final VCFFileReader reader = new VCFFileReader(plinkGeneratedVCF, false)) {
            final Map<String, Integer> sampleNameToOffset = reader.getFileHeader().getSampleNameToOffset();
            for (final VariantContext vc : reader) {
                for (final Genotype g : vc.getGenotypes()) {
                    Assert.assertEquals(g.isNoCall(), requestedMissingGenotypeIndices.contains(sampleNameToOffset.get(g.getSampleName())));
                }
            }
        }    
    }

}
