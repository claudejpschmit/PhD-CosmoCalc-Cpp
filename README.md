PhD-CosmoCalc-Cpp
=================

This is the C++ version of the PhD-CosmoCalc program, used to calculate basic cosmological functions and Fisher matrix analysis of the 21cm Brightness Temperature fluctuations. It utilizes the CAMB code, available at http://camb.info/. All necessary files are included in this project (due to compatibility requirements), so no external download is required. The code also requires ARES (http://ares.readthedocs.org/en/latest/index.html) for the 21cm fluctuations. Simply follow the installation guide on their website.

## Dependencies

The following libraries are used and should be installed before this code can be used:

	Armadillo C++ library,
	Boost,
	GSL_19,
	ARES (Accelerated Reionization Era Simulations code).


## Programming Style

We are generally following the [kernel style guide](https://www.kernel.org/doc/Documentation/CodingStyle), with the exception of indenting by 4, **not** 8, spaces.

## Building

A makefile for the application is provided.

### Building the application

In order to build the code, simply cd to the folder containing the Makefile and type:

    make calc

in order to generate the cosmology calculator with the excecutable 'calc'. 
To compile the analysis code that produces error ellipses from Fisher calculations, type:

    make analyse
    
The analyse program takes an integer that identifies the run for which the error ellipses shall be produces. Eg. one can run:

	./analyse 2

Simply typing:
    
    make

produces all excecutables.

### Building the documentation

EDIT: DOCUMENTATION BUILDING NOT UPDATED YET!

/To do that you need to be inside the build directory and make doc:
/
/    cd build
/    make doc
/

## Documentation

Using Doxygen. So write the documentation using [markdown](http://daringfireball.net/projects/markdown/syntax), which will enable easy HTML code generation. You can also find it helpful to [look at the Doxygen specific markdown docs](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html).


