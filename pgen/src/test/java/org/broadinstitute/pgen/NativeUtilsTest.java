/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

package org.broadinstitute.pgen;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

public class NativeUtilsTest {
    @Test
    public void testCreateTempDir() throws IOException {
        final String aPrefix = "aPrefix";
        final File tempDir = NativeLibraryUtils.createTempDir(aPrefix);
        Assert.assertTrue(tempDir.isDirectory());
        Assert.assertTrue(tempDir.getName().startsWith(aPrefix));
    }
  
}
