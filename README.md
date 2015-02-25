[![Build Status](http://apfel.mi.infn.it:9000/api/badge/github.com/scarrazza/apfel/status.svg?branch=master)](http://apfel.mi.infn.it:9000/github.com/scarrazza/apfel)

![alt text](https://github.com/scarrazza/apfel/raw/master/resources/logoapfel.png "Logo APFEL")

# APFEL: A PDF Evolution Library

Visit: http://apfel.hepforge.org and http://apfel.mi.infn.it/
 
APFEL is able to perform DGLAP evolution up to NNLO in QCD and to LO
in QED, both with pole and MSbar masses. The coupled DGLAP QCD+QED
evolution equations are solved in x-space by means of higher order
interpolations and Runge-Kutta techniques, allowing to explore
different options for the treatment of subleading terms.

The APFEL library implements the following tools:

- APFEL Web: APFEL provides an online web-app which integrates several
HEP softwares providing a complete suite for PDF analysis
(http://apfel.mi.infn.it/).

- APFEL GUI: APFEL provides the user with a graphical interface which
performs PDF, luminosity and DIS observables plots in real time,
granting an easy access to all the APFEL's functionalities.

## Download

You can obtain APFEL directly from the github repository:

https://github.com/scarrazza/apfel/releases

For the last development version you can clone the master code:

```Shell
git clone https://github.com/scarrazza/apfel.git
```

For the latest tag:

```Shell
git tag -l
git checkout tags/tag_name
```

## Installation 

Checkout the code and compile the code using the
following procedure:

```Shell
cd apfel
./configure --prefix=/where/install/apfel #(optional)
make && make install
```

By the default, if prefix is not set the program is installed in
/usr/local. If you define a custom prefix, remember to export
apfel/lib to the LD_LIBRARY_PATH. APFEL GUI requires ROOT (> 5), qmake
(> 4.5) and APFEL (> 1.0.1). The installation steps are:

```Shell
cd apfel/apfelGUI
qmake
make && make install
apfelgui
```
If you observe problems when running the GUI please set your LC_ALL=C.

## Contact Information

Maintainers: Valerio Bertone, Stefano Carrazza

Homepage: http://apfel.hepforge.org/
