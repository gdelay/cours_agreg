import numpy as np
import matplotlib.pyplot as plt



####### Figure 2 : on retrouve la fig 2 du sujet


# T : temps final
# N : nbre de pas de discretisation
# mu : parametre du modele
# Q : debit
# x0 : cond initiale
# zeta : cond initiale
def Euler_module(T, N, mu, Q, x0, zeta):
    X = np.zeros(N+1)
    Y = np.zeros(N+1)

    # conditions initiales
    X[0] = x0
    Y[0] = zeta

    # pas de temps
    h = T/N

    # methode d'Euler
    for i in range(N):
        X[i+1] = X[i] + h * (mu * X[i] * Y[i] - Q * X[i])
        Y[i+1] = Y[i] + h * (- mu * X[i] * Y[i] + Q * (zeta - Y[i]) )

    return [X,Y]


T = 20
dt = 0.02
N = 1000
mu = 1
zeta = 1
Qc = mu * zeta
x0 = 0.1
[X1,Y1] = Euler_module(T,N,mu, 0.8*Qc, x0, zeta)
[X2,Y2] = Euler_module(T,N,mu, 0.5*Qc, x0, zeta)
[X3,Y3] = Euler_module(T,N,mu, 0.3*Qc, x0, zeta)


tps = np.linspace(0,T,N+1)

plt.plot(tps,X1,'--')
plt.plot(tps,Y1,'-')
plt.text(14,0.82, 'Q/Qc = 0.8')
plt.text(14,0.15, 'Q/Qc = 0.8')

plt.plot(tps,X2,'--')
plt.plot(tps,Y2,'-')
plt.text(14,0.52, 'Q/Qc = 0.5')


plt.plot(tps,X3,'--')
plt.plot(tps,Y3,'-')
plt.text(14,0.72, 'Q/Qc = 0.3')
plt.text(14,0.25, 'Q/Qc = 0.3')

plt.ylabel('Concentrations')
plt.xlabel('Temps')

plt.show()


#######  Figure 3

# T : temps final
# N : nbre de pas de discretisation
# mu : parametre du modele
# Q : debit
# eps : epsilon
# z0 : cond initiale
def Euler_lac(T, N, mu, Q, eps, z0):
    Z = np.zeros(N+1)

    # conditions initiales
    Z[0] = z0

    # pas de temps
    h = T/N

    # methode d'Euler
    for i in range(N):
        Z[i+1] = Z[i] + h * Q * eps * (Q / mu - Z[i])

    return Z


T = 5000
N = 5000
eps = 0.01
mu = 1
z0 = 1
yobj = 0.5

Z1 = Euler_lac(T,N,mu,0.02,eps,z0)
Z2 = Euler_lac(T,N,mu,0.04,eps,z0)
Z3 = Euler_lac(T,N,mu,0.25,eps,z0)

Zconst = yobj * np.ones(N+1)


tps = np.linspace(0,T,N+1)

plt.plot(tps,Z1,'-', label='Q=0.02')
plt.plot(tps,Z2,'--', label='Q=0.04')
plt.plot(tps,Z3,':', label='Q=0.25')
plt.plot(tps,Zconst,'_',label='Concentration reglementraire')

plt.legend()
plt.show()


####### Figure 4

# T : temps final
# N : nbre de pas de discretisation
# mu : parametre du modele
# eps : epsilon
# z0 : cond initiale
def Euler_controle_lac(T, N, mu, eps, z0):
    Z = np.zeros(N+1)

    # conditions initiales
    Z[0] = z0

    # pas de temps
    h = T/N

    # methode d'Euler
    for i in range(N):
        Q = Z[i]/2
        Z[i+1] = Z[i] + h * Q * eps * (Q / mu - Z[i])

    return Z

T = 5000
N = 5000
eps = 0.01
mu = 1
z0 = 1
yobj = 0.1

Z1 = Euler_lac(T,N,mu,yobj/2,eps,z0)
Z2 = Euler_lac(T,N,mu,yobj/4,eps,z0)
Z3 = Euler_lac(T,N,mu,yobj,eps,z0)
Z4 = Euler_lac(T,N,mu,3*yobj,eps,z0)

Ztps = Euler_controle_lac(T, N, mu, eps, z0)

Zconst = yobj * np.ones(N+1)


tps = np.linspace(0,T,N+1)

plt.plot(tps,Z1,'-', label='Q=0.05')
plt.plot(tps,Z2,'--', label='Q=0.025')
plt.plot(tps,Z3,':', label='Q=0.1')
plt.plot(tps,Z4,'-.', label='Q=0.3')
plt.plot(tps,Ztps,'-.', label='Q=Q^*(t)')
plt.plot(tps,Zconst,'_',label='Concentration reglementraire')

plt.legend()
plt.show()
