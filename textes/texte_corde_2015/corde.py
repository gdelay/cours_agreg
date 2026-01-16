import numpy as np
import matplotlib.pyplot as plt

lambda_ = 30


def p(x):
    return -lambda_*np.exp(-(lambda_*x)**2)


def pder(x):
    return lambda_**3 * 2 * x * np.exp(-(lambda_*x)**2)


def k(x):
    return 1 + 0.8*np.sin(10*x) * np.cos(5*x)


def f(a, x):
    return p(x-a)


N = 100  # pas de discretisation

# affichage de la raideur de corde
X = np.linspace(0, 1, N)
K = [k(x) for x in X]
plt.plot(X, K)
plt.title("La tension k(x)")
plt.xlabel("x")
plt.ylabel("k(x)")
plt.show()


def MatM(N_):
    h = 1/(N_+1)
    ret = np.zeros((N_, N_))
    for i in range(1, N_):
        ret[i-1, i] = - k((i+0.5)*h) / h**2
    for i in range(1, N_+1):
        ret[i-1, i-1] = (k((i+0.5)*h) + k((i-0.5)*h)) / h**2
    for i in range(1, N_):
        ret[i, i-1] = -k((i+0.5)*h) / h**2

    return ret


def vectF(a, N_):
    h = 1/(N_+1)
    ret = np.zeros(N_)
    for i in range(0, N_):
        ret[i] = f(a, i*h)
    return ret


def solveY(a, N_):
    M_ = MatM(N_)
    F_ = vectF(a, N_)
    Y_ = np.linalg.solve(M_, F_)
    return Y_


ap = 0.666
# affichage de la charge pour a = ap
F = vectF(ap, N)
plt.plot(X, F)
plt.title("Charge pour a = 0.5")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()


# resolution et affichage pour a = ap
Y = solveY(ap, N)
plt.plot(X, Y)
plt.title("Deplacement pour a = 0.5")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.show()


# ~~~~~~~ Energie de rupture et sa derivee ~~~~~~~~~


def RN(a, N_):
    h = 1/(N_+1)
    Y_ = solveY(a, N_)

    E0 = 0.5*(k(0.5*h) * Y_[0] / h)**2
    EN = 0.5*(k((N_+0.5)*h) * Y_[N_-1] / h)**2
    return E0 + EN


def RNder(a, N_):
    h = 1/(N_+1)
    G = np.zeros(N_)
    for i in range(0, N_):
        G[i] = - pder((i+1)*h - a)

    M_ = MatM(N_)
    V = np.linalg.solve(M_, G)
    Y_ = solveY(a, N_)

    strength0 = k(0.5*h)**2 * Y_[0] * V[0] / h**2
    strengthN = k((N_+0.5)*h)**2 * Y_[N_-1] * V[N_-1] / h**2
    return strength0 + strengthN


# affichage de l'energie en fonction de a
M = 20
A = np.linspace(0.1, 0.9, M)
R = np.zeros(M)
for j in range(M):
    aj = A[j]
    R[j] = RN(aj, N)

plt.plot(A, R)
plt.show()

# ~~~~~~~ methode du grandient ~~~~~~~~~


def Gradient(N_, a0, rho, tol, nstop):
    a = a0
    n = 0
    while((abs(RNder(a, N_)) > tol) and n < nstop):
        a = a - rho * RNder(a, N_)
        n = n+1

    print("a = ", a, "n = ", n, "RNder = ", RNder(a, N_))
    return a


# ~~~~~~~ methode de Newton ~~~~~~~~~
# la methode de Newton semble ne pas marcher sur ce probleme


def Newton(N_, a0, tol, nstop):
    a = a0
    n = 0
    while((abs(RNder(a, N_)) > tol) and n < nstop):
        a = a - RN(a, N_) / RNder(a, N_)
        n = n+1

    print("a = ", a, "n = ", n, "RNder = ", RNder(a, N_))
    return a


A0 = 0.5
RHO = 0.008
TOL = 1e-6
NSTOP = 200

Gradient(N, A0, RHO, TOL, NSTOP)
Newton(N, A0, TOL, NSTOP)
