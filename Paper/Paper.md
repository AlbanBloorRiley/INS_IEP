---
title: 'INS_IEP: A Matlab package for fitting Inelastic Neutron Scattering data.'
tags:
  - Matlab 
  - Inelastic Neutron Scattering
  - Inverse Eigenvalue Problems
  - Deflation
  - NonLinear Least Squares
authors:
  - name: Alban Bloor Riley
    orcid: 0000-0001-9431-709X
    affiliation: 1
    corresponding: true
  - name: Marcus Webb
    orcid: 0000-0002-9440-5361
    affiliation: 1
  - name: Michael L. Baker
    orcid: 0000-0002-8246-3177
    affiliation: 2
affiliations:
  - name: The Uinversity of Manchester, Department of Mathematics
    index: 1
  - name: The Uinversity of Manchester, Department of Chemistry
    index: 2
date: June 2025
bibliography: Paper.bib

---

# Summary

`INS_IEP` is a MATLAB package used to fit experimental  inelastic neutron scattering (INS) data to a spin Hamiltonian model. INS experiments yield information about the eigenvalues of the spin Hamiltonian matrix that describes the sample being investigated, this can then be formulated as an associated inverse eigenvalue problem (IEP). `INS_IEP` uses deflated numerical optimisation techniques to find multiple spin systems that minimise this problem. 



# Statement of need


Inelastic neutron scattering (INS) is a spectroscopic technique used to measure the magnetic excitations in materials with interacting electron spins, such as single ions or molecular-based magnets. The  experiments are able to measure the energy between quantum spin states, these energy are then associated with the eigenvalues of the Hamiltonian matrix that describes the quantum spin dynamics of the compound in question [furrer_magnetic_2013; baker_neutron_2012; baker_spectroscopy_2014]. This information is crucial in understanding quantum phenomena and potentially can help  utilise electronic quantum spins in new quantum applications such as information sensing and processing. The particular parametrised inverse eigenvalue problem that arises from an INS experiment is determined by the Hamiltonian model that is used to describe the spin system, and the eigenvalues that have been found experimentally. This package is heavily inspired by the MATLAB packages `easyspin` [@stoll_easyspin_2006] for simulating and fitting Electron Paramagnetic Resonance (EPR) spectra, and `mint' [@Mint] for simulating INS spectra. `INS_IEP` inherits the  "Spin System" syntax from `easyspin` and `mint', meaning that they are all cross compatible. 



# Inverse Eigenvalue Problem 

The INS experiements provide eigenvalues of the Spin Hamiltonian matrix of the sample, the task of calculating the matrix from the eigenvalues is an inverse eigenvalue problem

Let $A(x)$ be the affine family of matrices,
$$A( x) = A_0 + \sum^\ell_{i=1} x_i A_i,$$
where $x\in\mathbb R^\ell$ and $A_0,\dots,A_\ell \in \mathbb R^{n\times n}$ are linearly independent symmetric matrices, and denote the ordered eigenvalues of $A(x)$ as $\lambda_1(x)\leq\dots\leq\lambda_n(x)$.
Then the least squares inverse eigenvalue problem is to find the parameters $x \in \mathbb R^\ell$ that minimises
$$  F(x,\rho) = \frac 1 2 ||r(x,\rho)||^2_2 = \frac 1 2 \sum^m_{i=1}(\lambda_{\rho_i}( x) - \lambda_i^*)^2 $$ , where  $\lambda_1^*\leq\ldots\leq\lambda_m^*$ are the experimental eigenvalues. 


Given real numbers $\lambda_1^*\leq\ldots\leq\lambda_m^*$, where $m \leq n$, find the parameters $x \in \mathbb R^\ell$ and the permutation $\rho$ $\in S_n$


# Overview of Methods

As far as we are aware this is the first time that the fitting of INS data has been explicitly formulated as an IEP. The number of eigenvalues that can be probed via INS experiments varies  depending on the equipment and sample in question, meaning that the fitting problem is often under (or even over) determined. The IEP is also highly nonlinear and due to the experimental nature of the data there is no guarentee that the problem is not ill-posed. One consequence of this is that the solution space my be very 'bumpy', that is there may exist many local minimisers to the problem. We seek to solve this problem by the use of Deflation, a numerical technique used to find multiple solutions to systems of equations [@farrell_deflation_2015; farrell_deflation_2020]. 
`INS_IEP` contains  a deflated Newton method, and two deflated Gauss-Newton methods [@DeflationPaper] for finding multiple minimising systems.

`INS_IEP` also contains a Riemannian Gradient descent method [@RGDLP Paper], inspired by the Lift and Projection method [@chen_least_1996], specifically designed for solving IEPs. This method is guaranteed to converge to a minimum, but only linearly, it is also currently not possible to apply deflation to this method.





# Acknowledgements

ABR thanks the University of Manchester for a Dean’s Doctoral Scholarship. MW thanks the Polish National Science Centre (SONATA-BIS-9), project no. 2019/34/E/ST1/00390, for the funding that supported some of this research. 

# References
