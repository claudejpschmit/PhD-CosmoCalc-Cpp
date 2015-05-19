PhD-CosmoCalc-Cpp
=================

This is the C++ version of the PhD-CosmoCalc program, used to calculate basic cosmological functions and Fisher matrix analysis of the 21cm Brightness Temperature fluctuations. It utilizes Lesgourgues et al's CLASS code, available at class-code.net. All necessary files are included in this project (due to compatibility requirements), so no external download is required. 

## Programming Style

We are generally following the [kernel style guide](https://www.kernel.org/doc/Documentation/CodingStyle), with the exception of indenting by 4, **not** 8, spaces.

## Building

A makefile for the application is provided.

### Building the application

In order to build the code, simply cd to the folder containing the Makefile and type:

    make calc

in order to generate the cosmology calculator with the excecutable 'calc'. To test whether CLASS works as intended, type:
    
    make class_test

This generates an excexutable that runs the standard CLASS excecutable as given from the code directly downloaded from class-code.net. To run the test, type:

    ./class_test lcdm.ini

This should not generate any errors.
Simply typing:
    
    make

produces both excecutables.

### Building the documentation

EDIT: DOCUMENTATION BUILDING NOT UPDATED YET!

/To do that you need to be inside the build directory and make doc:
/
/    cd build
/    make doc
/

## Documentation

Using Doxygen. So write the documentation using [markdown](http://daringfireball.net/projects/markdown/syntax), which will enable easy HTML code generation. You can also find it helpful to [look at the Doxygen specific markdown docs](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html).


