# pgen-jni

A Java PGEN writer, for writing [HTSJDK](https://github.com/samtools/htsjdk)
[VariantContext](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/writer/VariantContextWriter.java)
objects to a [plink2](https://www.cog-genomics.org/plink/2.0) [PGEN](https://www.cog-genomics.org/plink/2.0/input#pgen)
file. The Java writer component implements the [HTSJDK](https://github.com/samtools/htsjdk)
[VariantContext](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/variant/variantcontext/writer/VariantContextWriter.java)
interface. See the plink2 pgen [spec](https://github.com/chrchang/plink-ng/tree/master/pgen_spec) for more information about the PGEN file format.

The Java implementation uses an underlying native component, which is built here using a combination of local source files, plus some source
files taken directly from the plink-ng repo. NOTE: This project is based on, and uses (an alpha version of) that plink-ng pgen writing source.
See the notes below on how to update this project to use a newer version of plink-ng source.

## Supported runtime platforms:
- linux intel x64
- macos intel x64

MacOS aarm64 is not yet supported (limitation of the dev-nokee gradle native code/JNI plugin).

The generated artifacts are:
- pgen.jar (Java jar), which includes the native components:
  - libpgen.dylib (for macOS 64-bit Intel)
  - libpgen.so  (for Linux 64-bit Intel)

Since the inclusion of the native components in the **pgen-jni** jar is partly a manual process, see below for instructions on how to
publish the jar.
## Project Structure

There are two (gradle) sub-projects:
1. **pgen**: The Java/JNI layer, which consists of Java code, plus a very thin C++ implementation of the JNI methods. The C++ code does minimal
work, limited only to converting to/from JNI types, and handling/propagating exceptions. Otherwise this code delegates to code in the **pgen-lib**
sub-project (see below). This sub-project is also where the most extensive testing and validation takes place, via TestNG tests. The Java tests
require a plink2 executeable on the local device; this is used for various validation and concordance tests that are run as part of the test suite.
2. **pgen-lib** (not to be confused with the plink-ng **pgenlib** library, or the **pgenlib** namespace): A pure C/C++ layer that provides a convenient
C/C++ interface for use by the JNI layer in the **pgen** subproject, with C++ unit tests via boost. This layer has no Java/JNI dependencies, and
can be compiled, tested, and even deployed independently of the pgen Java components as C API over pgen writing. This layer is in turn composed of
two logical parts:
    1. the C/C++ layer that provides an API to the Java JNI code, and which is a port of the plink-ng python pgenlib writing code (specifically, the
    pgen writing code in `plink-ng/2.0/Python/src/pgenlib/pgenlib.pyx`). This layer in turn delegates to the underlying plink-ng **pgenlib** pgen
    writing C/C++ code.
    2. the **pgenlib** C/C++ code files that are copied (unmodified) from the plink-ng repo
  
  There are two C++ namespaces exposed to callers of the C++ code in **pgen-lib**:
  - **pgenlib** - the C/C++ callable types and functions that are used by the JNI layer (and implemented in the **pgen-lib** subproject)
  - **plink2** - the C/C++ types and functions that are part of the plink2 implementation, that are used by the pgenlib implementation

## Project Development
For development and test execution, it is recommended to treat these as two separate projects, using the `VSCode` editor for the Java
pgen project, and `CLion` or some other C++ compatible IDE for the **pgen-lib** project and boost tests.

The [dev-nokee](https://docs.nokee.dev/manual/jni-library-plugin.html) gradle plugin is used to build the entire project as a single unit
(currently this only builds the native components for the architecture on which is its running; either macos or Linux), and to run the
aggregate test suite (`./gradle clean test` will build and run both sets for tests for the current architecture).

## Building pgen-jni

Building requires a Java 17+ JDK, a C++11-compatible compiler, [boost](https://www.boost.org/doc/libs/1_80_0/libs/test/doc/html/index.html),
and a recent [plink-ng](https://www.cog-genomics.org/plink/2.0/) executeable (must be on the path). The **pgen-lib** C++ unit tests require boost
to be installed in `/usr/local/boost`.

### Docker File
Currently, the project only builds components for the architecture on which the build is running. A [Docker file]() is included in the repo to
make it easy to build and test the entire project, including the native component on Linux, and to publish the resulting artfacts.

## Updating to newer versions of plink-ng
The **pgen-lib** project uses the plink-ng **pgenlib** source directly. The git commit hash of the plink-ng code currently in use is stored in a file
called README_plink_commit_hash.md.

To update this project to a newer version of plink-ng source, do the following:
1. Pull the [plink-ng](https://www.cog-genomics.org/plink) repo into a folder that is parallel to the folder for this repo. Make sure the
master branch of plink-ng is checked out.
2. From the root of this project, run the `scripts\repoDiffs.sh` script. This will report any differences between the plink-ng code in this
repo, and the newer versions of those same files in the plink-ng repo. This step is advisory only - it provides a way to scan for any changes
that might be problematic.
3. Once you'satisfied that those changes are acceptable, run the `scripts\updatePlinkeCode.sh`. This will copy the new code from the plink-ng
project over the correspdonding files in this project.
4. The important part: Any changes that have been made to the python code in the plink-ng repo (as compared to the version of the python code
from the version of plink2 currently being used) need to be analyzed, and wherever appropriate, manually propagated (ported) to the C++ code in
**pgen-lib** (recall, the **pgen-lib** code is a port of the plink-ng python code, so if the python code changes, the ported C++ code may need
to change to reflect the changes; this is especially true if the plink-ng C++ code that is copied from the plink-ng pgenlib code has changed;
often these two layers change together.)

## Publishing/Releasing pgen-jni

The **pgen-jni** jar file that is published should include both a .dylib (for Macos) and a .so (for linux). The correct version is dynamically
selected at runtime based on the platform on which the jar file is running.

Because no attempt has been made in this project to support/enable cross-compilation of these components on a non-native platform, each
component must be built in it's own native environment. However, since we want both components included in the published jar, publishing
is only supported on linux (it is recommended to run the publishing task from a Docker container built from the Dockerfile included in
the repo), and the publishing protocol requires that the mac component be built first on Macos; with the resulting .dylib then checked
into the repo. The publish process can(/must) then be run on linux, where the checked-in mac component will be included in the published
jar.

The steps to publish are:
1. Build the mac component on macos (`./gradlew clean test jar`). Make sure all of the tests pass.
2. Copy the newly created `libpgen.dylib` file from `libs/main/macos` into the `mac_dylib` folder, overwriting the previously checked in version.
3. Commit the updated .dylib to the repo (create a single commit with only that change, in order to cleanly preserve the record of updates).
4. Push the changes to the repo.
5. Create a build environment on linux (it is easiest to use the Dockerfile in the root of the repo as the publish environment, since that
contains all of the dependencies necesssary to build the native linux component, the jar file, and to run the tests). If using the Docker
image, you mut also configure a git email and user name (`config --global user.email YOUR_EMAIL && git config --global user.name "YOUR NAME"`),
and the Artifactory user name and password must be set in the environment: `export ARTIFACTORY_USERNAME=username && export ARTIFACTORY_PASSWORD=password`).
6. Clone the repo and build the linux component (`./gradlew clean test jar`). Make sure all of the tests pass.
7. Create a tag for the release (`git tag your_tag -m "Your message."`)
8. Run the release/publish task (`./gradlew -Drelease=true clean publishAllPublicationsToArtifactoryRepository`). `-Drelease=true` is important
as it is the signal to the build to include the native mac component in the jar along with the native linux component.
9. Push the new tag up to the repo.

## Licensing
With the exception of the files in **pgen-lib** folder, this project is
licensed under the Apache License, Version 2.0 (the "License");
- You may not use this software except in compliance with the License.
- You may obtain a copy of the License at:

       http://www.apache.org/licenses/LICENSE-2.0

The files in the subproject located in the **pgen-lib** folder, some of which are from the plink-ng/plink2 project
(https://github.com/chrchang/plink-ng), are licensed under the Lesser GPL license. See the file **pgen-lib**/COPYING.lesser.
