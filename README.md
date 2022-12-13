# pgen-jni
Experimental jni wrapper of PGEN libraries

Work in Progress to provide PGEN Writer/Reader in java.

There are two (gradle) sub-projects:
  - pgen: the Java/JNI layer, which consists of Java code plus a thin C++ layer that delegates to pgen-lib,
	  with Java unit tests
  - pgen-lib: a pure C++ layer that is independent of JNI, with C++/boost unit tests 

The C++ unit tests require boost, which must be installed in /usr/local/boost.

