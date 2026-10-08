import numpy as np

a = 6*25.4
ell = 25.4*18
n = 1
R = 6.78/2
#R = 9.08/2

N = int(1e5)
x = 2*R*np.random.random(N) - R*np.random.random(N)
z = 2*R*np.random.random(N) - R*np.random.random(N)

mask = (x**2 + z**2 <= R**2)

x = x[mask]
z = z[mask]

I = np.cos(n*np.pi*z/ell)**2*np.cos(np.pi*x/a)**2
Aeff = np.mean(I)*np.pi*R**2
print(Aeff)