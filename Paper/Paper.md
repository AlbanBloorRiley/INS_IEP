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

Inelastic neutron scattering (INS) is a spectroscopic technique used to measure the magnetic excitations in materials with interacting electron spins. Fitting the experimental data to a spin Hamiltonian model can be formulated as an inverse eigenvalue problem (IEP). `INS_IEP` is a MATLAB package that uses deflated numerical optimisation techniques to find multiple solutions to this problem. The package requires and is fully compatible with `easyspin` [@stoll_easyspin_2006], a package for solving similar problems in electron paramagnetic resonance.

# Statement of need


Inelastic neutron scattering (INS) is a spectroscopic technique used to measure the magnetic excitations in materials with interacting electron spins, such as single ions or molecular-based magnets. The  experiments are able to measure the energy between quantum spin states, these energy are then associated with the eigenvalues of the Hamiltonian matrix that describes the quantum spin dynamics of the compound in question [@furrer_magnetic_2013; @baker_neutron_2012; @baker_spectroscopy_2014]. This information is crucial in understanding quantum phenomena and potentially can help  utilise electronic quantum spins in new quantum applications such as information sensing and processing. The particular parametrised inverse eigenvalue problem that arises from an INS experiment is determined by the Hamiltonian model that is used to describe the spin system, and the eigenvalues that have been found experimentally.



# The Spin Hamiltonian

The Spin Hamiltonian, $H$, is an approximation of the Hamiltonian that uses spin coordinates instead of orbital coordinates, it is  widely used to model data arising from many spectroscopy techniques @launay_electrons_2014]. It can be modeled as a linear combination of interaction terms, we will focus on the zero field interaction, $H_{ZFI}$, and the elctron-electron interaction, $H_{EEI}$:

$$H = H_{ZFI} + H_{EEI} $$

Both of these terms can themselves be modelled as the linear sum of other basis matrices. The zero field interaction can be written as:

$$H_{ZFI} = \sum_{-k\leq q \leq k} B^q_kO^q_k$$

where the $O^q_k$ are Stevens Operators [@], and $B^q_k$ the assosiated parameter. When there are multiple spin centres it is necessary to take kronecker products of the operator with identity matrices of the appropriate for each other spin centre.

When there are multiple spin centres it is also necessary to include an electron-electron interaction term, H_{EEI}. This term will be the sum of interaction terms between each pair of spin centres:

$$H_{EEI} = -\sum_{i\neq j} J_{ij}  S_i\cdot  S_j$$

where $ S_i$ is the vector of spin operators $ S_i = [S_x, S_y, S_z]$ for the $i$th spin centre, and $J_{ij}$ is the unknown parameter. Note that in the isotropic case $J$ can be thought of as a scalar value, but in the anisotripc case will be a matrix where the off diagonals are skew symmetric. While the summation is in theory over all spin centre comninations, in practice many of these contributions will be negligible - often only the nearest neighbour interactions are modeled. 


# Mathematical Background 
## Inverse Eigenvalue Problem


The INS experiements provide eigenvalues of the Spin Hamiltonian matrix of the sample, the task of calculating the matrix from the eigenvalues is an inverse eigenvalue problem:

Let $A(x)$ be the affine family of matrices,

$$A( x) = A_0 + \sum^\ell_{i=1} x_i A_i,$$

where $x\in\mathbb R^\ell$ and $A_0,\dots,A_\ell \in \mathbb R^{n\times n}$ are linearly independent symmetric matrices, and denote the ordered eigenvalues of $A(x)$ as $\lambda_1(x)\leq\dots\leq\lambda_n(x)$.
Then the least squares inverse eigenvalue problem is to find the parameters $x \in \mathbb R^\ell$ that minimises

$$F(x) = \frac 1 2 ||r(x)||^2_2 = \frac 1 2 \sum^m_{i=1}(\lambda_i(x) - \lambda_i^*)^2$$

where  $\lambda_1^*\leq\ldots\leq\lambda_m^*$ are the experimental eigenvalues. In the case of INS fitting the $A_i$ basis matrices will be a combination of Stevens operators and electron-electron exchange terms. The IEP described above is formulated as an least squares problem, this is due to the fact that the number of eigenvalues that can be probed by INS experiments is often a small subset of the full spectrum

As far as we are aware this is the first time that the fitting of INS data has been explicitly formulated as an IEP. The advantage of this formulation is that there are explicit formulas for the derivatives of $r(x)$, the first derivative (Jacobian) is:

$$ J_r(x) = \begin{pmatrix}
        q_1(x)^TA_1q_1(x)&\dots &q_1(x)^TA_\ell q_1(x)\\
        \vdots&\ddots&\vdots\\
        q_m(x)^TA_1q_m(x)&\dots& q_m(x)^TA_\ell q_m(x)
    \end{pmatrix}.$$
    
And the second derivative (Hessian) is:
    
$$ (H_{r})_{ij}   = 2\sum^m_{k=1} (\lambda_k-\lambda_k^*) \sum^m_{\substack{t=1\\\lambda_t\neq\lambda_k}} \frac{(q_t^TA_iq_k)(q_t^TA_jq_k)}{\lambda_k-\lambda_t}.$$
## Methods

Access to analytical forms of the derivatives means it is not necessary to use the finite difference approximation that other approaches use, making the three optimisation  methods that ``INS_IEP`` uses faster and more accurate. All of the methods used are iterative schemes of the form $x^{k+1} = x^k +p^k$ where the step $p^k$ uniquely defines each algorithm:

-  Newton's method: $p^k = (J_r^TJ_r + H_rr){-1}J_r^Tr$
-  Gauss-Newton method: $p^k = (J_r^TJ_r)^{-1}J_r^Tr$
-  Lift and Projection Method: $p^k = B^{-1}J_r^Tr$

The Lift and Projection method is a a Riemannian Gradient descent method [@RGDLP Paper], inspired by the Lift and Projection method [@chen_least_1996], specifically designed for solving IEPs. The matrix $B$ is a Gram matrix formed from the frobenius inner products of the basis matrices: $B_{ij} = \langle A_i, A_j\rangle_F$. In [@RGDLP Paper] it is proven that the method is a strictly descending algorthim, that reduces the value of the objective function every step.

## Deflation 

The number of eigenvalues that can be probed via INS experiments varies  depending on the equipment and sample in question, meaning that the fitting problem is often under (or even over) determined. The IEP is also highly nonlinear and due to the experimental nature of the data there is no guarentee that the problem is not ill-posed. One consequence of this is that the solution space my be very 'bumpy', that is there may exist many local minimisers to the problem. For example in figure \autoref{fig:mn12}, there are clearly 4 distinct solutions (for more details see Example1_Mn12 in the examples folder). We seek to solve the problem of multiple local minima by the use of Deflation, a numerical technique used to find multiple solutions to systems of equations [@farrell_deflation_2015; farrell_deflation_2020]. Fortunately it is cheap to apply deflation for the above methods, it is simply a change to the length of the step - notably this means that the direction of each step does not change [@Deflation_Paper] . 

<!--![Left: Contour plot of how $F$ varies with the two parameters $B^2_2$ and $B_4^4$ for the molecule Mn\_12 as described in [@bircher_transverse_2004]. Right: Convergence plot for finding each minimum using deflation 
\label{fig:mn12}](Mn12_figure.png) 
-->
![Contour plot of how $F$ varies with the two parameters $B^2_2$ and $B_4^4$ for the molecule Mn\_12 as described in [@bircher_transverse_2004]\label{fig:contour}](Mn12_contour.png){ width=400 }
![Rate of converge for each deflation.\label{fig:convergence}](Mn12_convergence.png){ width=400 }

# Acknowledgements

ABR thanks the University of Manchester for a Dean’s Doctoral Scholarship. MW thanks the Polish National Science Centre (SONATA-BIS-9), project no. 2019/34/E/ST1/00390, for the funding that supported some of this research. 

# References


