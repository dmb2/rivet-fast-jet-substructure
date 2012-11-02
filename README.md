Hacked Rivet Fast Jet library 
=====================

This is a hacked version of the Fast-Jet library shipped with
[Rivet](rivet.hepforge.org). It has no meaning without having Rivet installed.

## Introduction

The BOOST 2012 working group has identified the following substructure
observables:
* Trimming/mass-drop/pruning/filtering
* Template method
* N-subjettiness
* Jet charge
* Angular correlation function
* Diploarity

Currently, Jet Charge and 2-subjettiness have been implemented. 

### Jet Charge
Jet Charge is defined as
```
$$
Q_{\kappa}^i = \frac{1}{(P_T^{jet})^\kappa}\sum_{j\in jet}Q_j(p_T^j)^\kappa
$$
```
Put in words, its the sum of the charge of the constituent particles
weighted by the particle's transverse momentum. 
### Jet Dipolarity
Jet dipolarity is a p_T weighted sum over the radii of sub*jets. It is
typically only defined for jets with two subjets. 

## Installing
Copy the folders include/ and src/ to 
/path/to/rivet/build/rivet/
```
cd /path/to/rivet/build/
make 
make install
```
Now profit from substructure and measure your favorite Higgs coupling.

_NB_ install_hacked_fastjets.sh will do the install/build process, but it is
very dumb so read it before you use it... no guarentees if it nukes your rivet
build.