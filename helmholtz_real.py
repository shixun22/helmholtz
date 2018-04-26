"""
Helmholtz-Hodge decomposition using Fourier transform
*** Not using np.fft, slower, but takes careful account of Nyquist frequency when N is even. => Guarantees real decomposed V fields, and zero divergence of V_solenoidal ***
Xun Shi 26.04.2018
--------------------------------------
Given 3D velocity field, decompose to solenoidal and compressive parts

Note:
1. Only for uniform grid with dx = dy = dz. Extensions are easy to implement.

2. Helmholtz-Hodge decomposition performed with spectral method should only apply to relatively smooth fields, i.e. with little power on small scales

3. For even NX, NY, NZ, decomposed fields can be complex, with the imaginary part coming from the real part of the kmode at Nyquist frequency. In principle, when doing first derivatives, the Nyquist frequency kmode should be dropped so as not to break symmetry. See e.g. footnote on page 4 of 

http://math.mit.edu/~stevenj/fft-deriv.pdf

However, when the field is smooth enough, the imaginary part caused by the Nyquist frequency kmode should be negligible.   
"""

import numpy as np
#import matplotlib.pyplot as plt

V = np.fromfile('v0_xyz_16_16_128.dat').reshape(3,16,16,128)

if V.shape[0]!=3: print 'Not a 3D field in required form!'

NX, NY, NZ = V.shape[1:]

def fftfreq_doublenyquist(N):
    """
    fftfreq_ has only one negative Nyquist frequency -0.5
    fftfreq_doublenyquist has both -0.5 and 0.5
    """
    return np.arange(-(N/2), (N/2)+1) * 1. / N

x = np.arange(NX)
y = np.arange(NY)
z = np.arange(NZ)
kx = fftfreq_doublenyquist(NX)
ky = fftfreq_doublenyquist(NY)
kz = fftfreq_doublenyquist(NZ)
k2 = kx[:,np.newaxis,np.newaxis]**2 + ky[np.newaxis,:,np.newaxis]**2 + kz[np.newaxis,np.newaxis,:]**2
k2[np.where(abs(k2)==0.)] = 1.

V_k = np.array([[[[(V[xyz] * np.exp(-2.* np.pi * 1j * (kx[ikx] * x[:,np.newaxis,np.newaxis] + ky[iky] * y[np.newaxis,:,np.newaxis] + kz[ikz] * z[np.newaxis,np.newaxis,:]))).sum() for ikz in range(kz.size)] for iky in range(ky.size)] for ikx in range(kx.size)] for xyz in range(3)])

if (NX%2)==0:
    V_k[:,0,:,:] = V_k[:,0,:,:]/2.
    V_k[:,-1,:,:] = V_k[:,-1,:,:]/2.

if (NY%2)==0:
    V_k[:,:,0,:] = V_k[:,:,0,:]/2.
    V_k[:,:,-1,:] = V_k[:,:,-1,:]/2.

if (NZ%2)==0:
    V_k[:,:,:,0] = V_k[:,:,:,0]/2.
    V_k[:,:,:,-1] = V_k[:,:,:,-1]/2.

## recontruct V from inverse Fourier transformation
#V_reconstruct = np.array([[[[(V_k[xyz] * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)] for xyz in range(3)])

div_V_k = (V_k[0] * kx[:,np.newaxis,np.newaxis] + V_k[1] * ky[np.newaxis,:,np.newaxis] + V_k[2] * kz[np.newaxis,np.newaxis,:]) / k2

V_solenoidal_k = np.array([V_k[0] - div_V_k * kx[:,np.newaxis,np.newaxis], V_k[1] - div_V_k * ky[np.newaxis,:,np.newaxis], V_k[2] - div_V_k * kz[np.newaxis,np.newaxis,:]])

V_solenoidal = np.array([[[[(V_solenoidal_k[xyz] * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)] for xyz in range(3)])


# check if the solenoidal part is real 
print 'max of imaginary part of V_solenoidal:', abs(V_solenoidal.imag).max()

# check if the solenoidal part really divergence-free
divVs_r = np.array([[[((V_solenoidal_k[0] * kx[:,np.newaxis,np.newaxis] + V_solenoidal_k[1] * ky[np.newaxis,:,np.newaxis] + V_solenoidal_k[2] * kz[np.newaxis,np.newaxis,:]) / k2 * 2. * np.pi * 1j * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)])

print 'max of divergence(V_solenoidal):', abs(divVs_r).max()

## plot one slice of the decomposed field on X-Y plane
#X, Y = np.meshgrid(range(NY), range(NX))
#plt.figure()
#plt.quiver(X, Y, V_solenoidal[0][:,:,0], V_solenoidal[1][:,:,0])
#plt.ion()
#plt.show()

