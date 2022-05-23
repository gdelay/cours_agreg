import numpy as np
import matplotlib.pyplot as plt

lambda_ = 30


def p(x):
    return -lambda_*np.exp(-(lambda_*x)**2)


def pder(x):
    return lambda_**3 * 2 * x * np.exp(-(lambda_*x)**2)


def k(x):
    return 1 + 0.8*np.sin(10*x) * np.cos(5*x)


def f(a,x):
    return p(x-a)


N = 100  # pas de discretisation

# affichage de la raideur de corde
X = np.linspace(0,1,N)
K = [k(x) for x in X]
plt.plot(X,K)
plt.title("La tension k(x)")
plt.xlabel("x")
plt.ylabel("k(x)")
plt.show()


def MatM(N):
    h = 1/(N+1)
    ret = np.zeros((N,N))
    for i in range(1,N):
        ret[i-1,i] = - k( (i+0.5)*h ) / h**2
    for i in range(1,N+1):
        ret[i-1,i-1] = (k( (i+0.5)*h ) + k( (i-0.5)*h ))/ h**2
    for i in range(1,N):
        ret[i,i-1] = -k( (i+0.5)*h )/ h**2

    return ret


def vectF(a,N):
    h = 1/(N+1)
    ret = np.zeros(N)
    for i in range(0,N):
        ret[i] = f(a,i*h)
    return ret


def solveY(a,N):
    M = MatM(N)
    F = vectF(a, N)
    Y = np.linalg.solve(M, F)
    return Y

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


####### Energie de rupture et sa derivee #######

def RN(a,N):
    h = 1/(N+1)
    Y = solveY(a, N)

    return 0.5*(k(0.5*h) * Y[0] / h)**2 + 0.5*(k((N+0.5)*h) * Y[N-1] / h)**2


def RNder(a,N):
    h = 1/(N+1)
    G = np.zeros(N)
    for i in range(0, N):
        G[i] = - pder((i+1)*h - a)

    M = MatM(N)
    V = np.linalg.solve(M, G)
    Y = solveY(a, N)

    return k(0.5*h)**2 * Y[0] * V[0] / h**2 + k((N+0.5)*h)**2 * Y[N-1] * V[N-1] / h**2


# affichage de l'energie en fonction de a
M = 20
A = np.linspace(0.1,0.9,M)
R = np.zeros(M)
for i in range(M):
    ai = A[i]
    R[i] = RN(ai, N)

plt.plot(A, R)
plt.show()

####### methode du grandient  #######



def Gradient(N, a0, rho, tol, nstop):
    a = a0
    n = 0
    while((abs(RNder(a, N)) > tol) and n < nstop):
        a = a - rho * RNder(a, N)
        n = n+1

    print("a = ", a, "n = ", n, "RNder = ", RNder(a, N))
    return a


####### methode de Newton  ######
##### la methode de Newton semble ne pas marcher sur ce probleme

def Newton(N, a0, tol, nstop):
    a = a0
    n = 0
    while((abs(RNder(a, N)) > tol) and n < nstop):
        a = a - RN(a, N) / RNder(a, N)
        n = n+1

    print("a = ", a, "n = ", n, "RNder = ", RNder(a, N))
    return a


a0 = 0.5
rho = 0.008
tol = 1e-6
nstop = 200

Gradient(N, a0, rho, tol, nstop)
Newton(N, a0, tol, nstop)
