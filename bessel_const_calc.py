import numpy as np
import matplotlib.pyplot as plt

Aeff = []

for n in [1, 3, 5, 7, 9]:

    N = 100000000

    R = 4.5

    x = 2*R*np.random.random(N) - R
    z = 2*R*np.random.random(N) - R

    r = np.sqrt(x**2 + z**2)

    a = 6*25.4
    l = 18*25.4


    r = r[r < R]


    epsr = 2.54

    k0 = np.sqrt((np.pi/a)**2 + (n*np.pi/l)**2)
    kperp = (np.pi/a)
    kperp_in = np.sqrt(kperp**2 + (k0**2)*(epsr - 1))

    def J0(x):
        J = 1 - (x/2)**2 +(1/4)*(x/2)**4 - (x/2)**6/(36)
        return J

    E0 = np.cos(np.pi*x[(x**2 + z**2) < R**2]/a)*np.cos(n*np.pi*z[(x**2 + z**2) < R**2]/l)


    I = (J0(kperp*R)/J0(kperp_in*R))*np.mean(E0*J0(kperp_in*r))*np.pi*R**2
    I2 = np.mean(E0**2)*np.pi*R**2
    Aeff.append(I)
    print(I, I2, I/I2)

plt.plot(Aeff)
plt.show()

