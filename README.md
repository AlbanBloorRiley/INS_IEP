# INS_IEP


## Installation
* Download the package on your PC.
* Dowload the 'EasySpin' package - details below.
* Open MATLAB
* Add both the 'INS_IEP' and 'easyspin' directories to your MATLAB path
* Go to the directory 'INS_IEP'


### Dependencies 
Note that this package requires [`easyspin`](https://easyspin.org) to be installed and added to the MTLAB path. This can be downladed at https://easyspin.org.

### Test
You are then ready to run the tests, by openning the Tests/ directory and then running `iitest`


## Repository Organization

The INS_IEP repository is organized as follows:

Component Formation/        source code for forming the various matrices and functions used by INS_IEP
Conversions/      source code for converting between easyspin sytax and mathmatical optimisation syntax
Examples/     collection of INS_IEP example scripts
Numerical Methods/    source code for the numerical methods used
Paper/    source code for the JOSS paper (under review)
Tests/         unit tests for INS_IEP
