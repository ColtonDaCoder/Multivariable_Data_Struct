import numpy as np

lam = 210*1e-9
c0 = 299792458
#n = (1)**0.5

omega = (2*np.pi*c0)/(lam)
eps0 = 8.854*1e-12
norm_A_sqr = 2
C0 = np.abs((eps0*omega)/(2))
print(C0)


