# pgen-jni

A Java PGEN writer, for writing [HTSJDK](https://github.com/samtools/htsjdk)
[VariantContext](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/writer/VariantContextWriter.java)
objects to a [plink2](https://www.cog-genomics.org/plink/2.0) [PGEN](https://www.cog-genomics.org/plink/2.0/input#pgen)
file. The Java writer component implements the [HTSJDK](https://github.com/samtools/htsjdk)
[VariantContext](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/writer/VariantContextWriter.java)
interface. See the plink2 pgen [spec](https://github.com/chrchang/plink-ng/tree/master/pgen_spec) for more information about the PGEN file format.

The Java implementation uses an underlying native component, which is built here using a combination of local
source files, plus some source files taken directly from the plink2 (that are used by plink2 to build the pgen-lib target).
## Artifacts

### Supported runtime platforms:
- linux x64
- macos x64

aarm64 is not yet supported (limitation of the dev-nokee gradle native code/JNI plugin)

The generated artifacts are:
- pgen.jar (Java jar)
  - libpgen.dylib (for macOS 64-bit)
  - libpgen.so  (for Linux 64-bit Intel1)
### Building pgen-jni

Building requires a Java 17+ JDK, a C++11-compatible compiler, [boost](https://www.boost.org/doc/libs/1_80_0/libs/test/doc/html/index.html),
and a recent [plink2](https://www.cog-genomics.org/plink/2.0/) executeable (must be on the path). The pgen-lib C++ unit tests require boost
to be installed in /usr/local/boost.

Currently, the project only builds components for the architecture on which the build is running. A [Docker file]() is provided to
make it easy to build the native component on Linux and publish the results.
### Publishing/Releasing pgen-jni

The pgen-jni jar file that is published should include both a .dylib (for Macos) and a .so (for linux). The correct version is dynamically
selected at ruimte based on the platform on which the jar file is running.

Because no attempt has been made in this project to support/enable cross-compilation of these components on anything other than their own
native environment(s), each component must be built in it's own native environment. However, since we want both components included in
the published jar, publishing is only supported on linux (it is recommended to run the publishing task from a Docker container built from
the Dockerfile included in the repo), and the publishing protocol requires that the mac component be built first on Macos; with the resulting
.dylib then checked into the repo. The publish process can then be run on linux.

The steps to publish are:
1. Build the mac component on macos (./gradlew clean test jar). Make sure the tests all succeed.
2. Copy the newly created libpgen.dylib file from libs/main/macos into the mac_dylib folder, overwriting the previously checked in version.
3. Commit the updated .dylib to the repo (create a single commit with only that change, in order to cleanly preserve the record of updates).
4. Push the changes to the repo.
5. Create a build environment on linux (it is easiest to use the Dockerfile in the root of the repo as the publish environment, since that
contains all of the dependencies necesssary to build the native linux component, the jar file, and to run the tests). If using the Docker
image, you mut also configure a git email and user name (config --global user.email YOUR_EMAIL && git config --global user.name "YOUR NAME"),
and the Artifactory user name and password must be set in the environment: ()"export ARTIFACTORY_USERNAME=username && export ARTIFACTORY_PASSWORD=password").
6. Clone the repo and build the linux component (./gradlew clean test jar). Make sure all of the tests have passed.
7. Create a tag for the release (git tag your_tag -m "Your message.")
8. Run the release/publish task (./gradlew -Drelease=true clean publishAllPublicationsToArtifactoryRepository)
9. Push the new tag up to the repo.

### Projects

There are two (gradle) sub-projects:
- pgen: The Java/JNI layer, which consists of Java code, plus a very thin C++ layer that implements the JNI methods. This code does minimal
work, limited to converting to/from JNI types and handling/propagating exeptions, but otherwise delegating to to pgen-lib. This is also
where the most extensive testing and validation takes place, via TestNG tests.
- pgen-lib: A pure C/C++ layer, with C++/boost unit tests. This layer has no Java/JNI 
dependencies, and can be compiled, tested, and deployed independently of the pgen Java or native components defined in PGEN.

The two sub-projects have separate tests (TestNG for the pgen Java layer, and C++/boost for the pgen-lib unit tests).
For development and test execution, it is recommended to treat these as two separate projects, using the VSCode editor for the Java
pgen project, and CLion or other C++ compatible IDE for the pgen-lib project and boost tests.

The [dev-nokee](https://docs.nokee.dev/manual/jni-library-plugin.html) gradle plugin is used to build the entire project as a single unit
(currently this only builds the native components for the architecture on which is its running; either macos or Linux), and to run the
aggregate test suite (`./gradle clean test` will build and run both sets for tests for the current architecture).

#### pgen
The Java tests require the presence of a plink2 executeable on the local device; this is used for varius validation and concordance
tests that are run as part of the test suite.

##### pgen-lib
There are two C++ namespaces exposed to callers of the C++ code in pgen-lib:
- pgenlib - the C/C++ calleable types and functions that are used by the JNI layer (and implemented in the pgen-lib subproject)
- plink2 - the C/C++ types and functions that are part of the plink2 implementation, that are used by the pgenlib implementation

## Licensing
With the exception of the files in pgen-lib folder, this project is
licensed under the Apache License, Version 2.0 (the "License");
- You may not use this software except in compliance with the License.
- You may obtain a copy of the License at:

       http://www.apache.org/licenses/LICENSE-2.0

The files in the subproject located in the pgen-lib folder, some of which are from the plink-ng/plink2 project
(https://github.com/chrchang/plink-ng), are licensed under the Lesser GPL license. See the file pgen-lib/COPYING.lesser.
