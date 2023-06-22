/**
 * Copyright (c) 2023, Broad Institute, Inc. All rights reserved.
 */

package org.broadinstitute.pgen;

import org.apache.commons.compress.utils.IOUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.SystemUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;

/**
 * Utilities to provide architecture-dependent native functions
 * 
 * IMPORTANT: be careful when modifiying this code, since most of it is used to extract a resource from a jar,
 * and as a result is neither covered by tests, nor ever executed by this project. It is only used to load
 * the native components when the jar produced by this project is used in another project.
 */
public final class NativeLibraryUtils {

    private NativeLibraryUtils(){}

    private record Resource(String path, Class<?> relativeClass) {
        /**
         * Get the contents of this resource as an InputStream
         * @throws IllegalArgumentException if resource cannot be read
         * @return an input stream that will read the contents of this resource
         */
        public InputStream getResourceContentsAsStream() {
            //final Class<?> clazz = getRelativeClass();

            final InputStream inputStream;
            if (relativeClass == null) {
                inputStream = ClassLoader.getSystemResourceAsStream(path);
                if (inputStream == null)
                    throw new IllegalArgumentException("Resource not found: " + path);
            } else {
                inputStream = relativeClass.getResourceAsStream(path);
                if (inputStream == null)
                    throw new IllegalArgumentException("Resource not found relative to " + relativeClass + ": " + path);

            }

            return inputStream;
        }

        /**
         * Writes the an embedded resource to a file.
         * File is not scheduled for deletion and must be cleaned up by the caller.
         * @param resource Embedded resource.
         * @param file File path to write.
         */
        @SuppressWarnings("deprecation")
        public void writeResource(File file) {
            InputStream inputStream = getResourceContentsAsStream();
            OutputStream outputStream = null;
            try {
                outputStream = FileUtils.openOutputStream(file);
                org.apache.commons.io.IOUtils.copy(inputStream, outputStream);
            } catch (IOException e) {
                throw new RuntimeException(String.format("Unable to copy resource '%s' to '%s'", path, file), e);
            } finally {
                org.apache.commons.io.IOUtils.closeQuietly(inputStream);
                org.apache.commons.io.IOUtils.closeQuietly(outputStream);
            }
        }
    }

    /**
     * @return true if we're running on a Mac operating system, otherwise false
     */
    public static boolean runningOnMac() {
        return SystemUtils.IS_OS_MAC;
    }

    /**
     * @return true if we're running on a Linux operating system, otherwise false
     */
    public static boolean runningOnLinux() {
        return SystemUtils.IS_OS_LINUX;
    }

    /**
     * Loads a native library stored on our classpath by extracting it to a temporary location and invoking
     * {@link System#load} on it
     *
     * @param libraryPathInJar absolute path to the library file on the classpath
     * @return true if the library was extracted and loaded successfully, otherwise false
     */
    public static boolean loadLibraryFromClasspath( final String libraryPathInJar ) {
    //    Utils.nonNull(libraryPathInJar);
    //    Utils.validateArg(libraryPathInJar.startsWith("/"), "library path in jar must be absolute");

        try {
            final File extractedLibrary = writeTempResourceFromPath(libraryPathInJar, NativeLibraryUtils.class);
            extractedLibrary.deleteOnExit();
            System.load(extractedLibrary.getAbsolutePath());
        }
        catch ( UnsatisfiedLinkError e ) {
            return false;
        }

        return true;
    }

    /**
     * Create a resource from a path and a relative class, and write it to a temporary file.
     * If the relative class is null then the system classloader will be used and the path must be absolute.
     * The temporary file is automatically scheduled for deletion on exit.
     * @param resourcePath Relative or absolute path to the class.
     * @param relativeClass Relative class to use as a class loader and for a relative package.
     * @return a temporary file containing the contents of the resource, which is automatically scheduled
     * for deletion on exit.
     */
    private static File writeTempResourceFromPath(final String resourcePath, final Class<?> relativeClass) {
        //Utils.nonNull(resourcePath, "A resource path must be provided");
        try {
            final Resource resource = new Resource(resourcePath, relativeClass);
            final File tmpDir = createTempDir("nativeResource");
            final File tempFile = File.createTempFile(FilenameUtils.getBaseName(resource.path()), FilenameUtils.getExtension(resource.path()), tmpDir);
            tempFile.deleteOnExit();

            resource.writeResource(tempFile);
            return tempFile;
        } catch (final IOException e) {
            throw new RuntimeException("failure to write resource", e);
        }
    }

    /**
     * Creates a temp directory with the given prefix.
     *
     * The directory and any contents will be automatically deleted at shutdown.
     *
     * This will not work if the temp dir is not representable as a File.
     *
     * @param prefix       Prefix for the directory name.
     * @return The created temporary directory.
     */
    //@VisibleForTesting
    static File createTempDir(String prefix) {
        try {
            final Path tmpDir = Files.createTempDirectory(prefix).normalize();
            tmpDir.toFile().deleteOnExit();
            return tmpDir.toFile();
        } catch (final IOException | SecurityException e) {
            throw new RuntimeException(String.format(
                "Bad tmp dir: %s", e.getMessage()), e);
        }
    }
}
