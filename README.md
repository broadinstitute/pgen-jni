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

### Supported platforms:
- linux 64
- macos 64

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
make it easy to build the native components on Linux.
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





