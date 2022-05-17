#################### IMPORTS ##############
import numpy as np
import matplotlib.pyplot as plt
##################  PARAMETERS  ###########

N = 40  # number of time  sub-divisions
M = 50  # number of space sub-divisions

N_int = N/2  # intermediate time (plot)

T = 0.4  # final time
c = 0.2  # wave celerity


# initial position
def u0(x):
    if(x <= 0.3 or x >= 0.7):
        return 0
    if(x < 0.5):
        return x-0.3

    return 0.7-x


# initial velocity
def u1(x):
    return 0


# source term
def f(t,x):
    return 0


##################  MAIN CODE  ############

# steps
ht = T/N
hx = 1./M

# X = (x_j)
X = np.linspace(0, 1, M+1)

################  SCHEMA CENTRE  ##########

## initialisation of the solution
U = np.zeros(M+1)
for j in range(M+1):
    U[j] = u0(X[j])

U_ini = np.zeros(M+1)
for j in range(M+1):
    U_ini[j] = u0(X[j])

U_prev = np.zeros(M+1)
for j in range(M+1):
    U_prev[j] = U[j]

# first time step
for j in range(M+1):
    U[j] = U_prev[j] + ht * u1(X[j])


## time loop
CFL = c * ht / hx

print("CFL = ", CFL)

U_prev_prev = np.zeros(M+1)
U_int = np.zeros(M+1)

for n in range(2, N):
    tn = n*ht
    for j in range(M+1):
        U_prev_prev[j] = U_prev[j]
        U_prev[j] = U[j]

    for j in range(1, M):
        diff_term = CFL * CFL * (U_prev[j+1] - 2*U_prev[j] + U_prev[j-1])
        U[j] = 2*U_prev[j] - U_prev_prev[j] + diff_term #+ ht*ht*f(tn, X[j])
    U[0] = 0  # Dirichlet boundary condition
    U[M] = 0  # Dirichlet boundary condition

    if(n == N_int):  # save intermediate result
        for j in range(M+1):
            U_int[j] = U[j]

# exact solution (if the wave does not touch the boundary and if f=0 and u1 = 0)
ex_sol = np.zeros(M+1)
for j in range(M+1):
    ex_sol[j] = (u0(X[j]+c*T) + u0(X[j]-c*T))/2

## plot final solution
plt.plot(X,U, label="t=0.4")
plt.plot(X,U_int, label="t=0.2")
plt.plot(X,ex_sol, label="solution exacte t=0.4")
# plt.plot(X,U_ini, label="t=0")
plt.title("Solution numerique de l'equation des ondes")
plt.xlabel("x")
plt.ylabel("solution")
plt.legend()
plt.show()
