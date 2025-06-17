#################### IMPORTS ##############
import numpy as np
import matplotlib.pyplot as plt
##################  PARAMETERS  ###########

T = 0.4  # final time
c = 1  # wave celerity

# exact solution
def sol(t,x):
    return np.sin(np.pi*x) * np.cos(c*np.pi*t)

# initial position
def u0(x):
    return np.sin(np.pi*x)


# initial velocity
def u1(x):
    return 0


# source term
def f(t,x):
    return 0


##################  Error computation  ############
def error(N,M):
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
        U[j] = U_prev[j] + ht * u1(X[j]) # passer ce terme a l'ordre superieur

    CFL = c * ht / hx
    
    print("CFL = ", CFL)
    
    U_prev_prev = np.zeros(M+1)
    U_int = np.zeros(M+1)
    
    ## main loop (over time)
    for n in range(2, N):
        tn = n*ht
        for j in range(M+1):
            U_prev_prev[j] = U_prev[j]
            U_prev[j] = U[j]
        
        for j in range(1, M):
            diff_term = CFL * CFL * (U_prev[j+1] - 2*U_prev[j] + U_prev[j-1])
            U[j] = 2*U_prev[j] - U_prev_prev[j] + diff_term + ht*ht*f(tn, X[j])
        U[0] = 0  # Dirichlet boundary condition
        U[M] = 0  # Dirichlet boundary condition

    ## exact solution (if the wave does not touch the boundary and if f=0 and u1 = 0)
    ex_sol_T = np.zeros(M+1)
    for j in range(M+1):
        ex_sol_T[j] = sol(T,X[j])

    ## error computation at final time
    err_T = np.zeros(M+1)
    for j in range(M+1):
        err_T[j] = np.abs(sol(T,X[j]) - U[j])
    err_max_T = np.max(err_T)
    return err_max_T


##################  Convergence rates  ############

## zone a deboguer

meshes = [25,50,100,200,400]
Err = []
liste_h = []
CFL = 0.5
for M in meshes:
    ## corresponding time discretization
    N = T*CFL/(c*M)
    print(N)
    
    ## error
    err = error(N,M)
    Err.append(err)
    liste_h.append(1./M)
    print(err)


## plot the convergence errors in log scale
Err_log = np.log10(Err)
h_log = np.log10(liste_h)
plt.plot(h_log, Err_log)
plt.xlabel("log(h)")
plt.ylabel("log(erreur)")
plt.title("Etude de convergence pour conditions de Dirichlet")
plt.show()


## polyfit
z = np.polyfit(h_log,Err_log,1)
print('z = ', z)
## on retrouve numeriquement que la methode est d'ordre 2
