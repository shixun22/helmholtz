# helmholtz
Helmholtz-Hodge decomposition using fft

Contains:
helmholtz.py  -- python code for the decomposition
v0_xyz_16_16_128.dat  -- sample vector field to be decomposed

--------------------------------------
Given 3D vector field, decompose to solenoidal and compressive parts

Note:
1. Only for uniform grid with dx = dy = dz. Extensions are easy to implement.

2. Helmholtz-Hodge decomposition performed with spectral method should only apply to relatively smooth fields, i.e. with little power on small scales

3. For even NX, NY, NZ, decomposed fields can be complex, with the imaginary part coming from the real part of the kmode at Nyquist frequency. In principle, when doing first derivatives, the Nyquist frequency kmode should be dropped so as not to break symmetry. See e.g. footnote on page 4 of 

http://math.mit.edu/~stevenj/fft-deriv.pdf

However, when the field is smooth enough, the imaginary part caused by the Nyquist frequency kmode should be negligible. 


