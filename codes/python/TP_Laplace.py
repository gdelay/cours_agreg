#################### IMPORTS ##############
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import numpy.linalg as LA
##################  PARAMETERS  ###########

M = 200 # number of space sub-divisions


## RHS function
def f(x):
    #return x
    return x*x


def f_mixte(x):
    return x*x


## Exact solution
def sol_ex(x):
    #return -x**3/6. + x/6.
    return x/12. * (1.-x*x*x) + x/12.


def sol_mixte(x):
    return x/3. * (1.-x*x*x/4.)


##################  MAIN CODE  ############

## Dirichlet data
alpha = sol_ex(0)
beta = sol_ex(1)
## note : on retrouve les conditions de Dirichlet homogenes en posant alpha = beta = 0



#############  DIRICHLET BC  ############

def Dirichlet(M,rhs_fun, A, B):
    ## step
    h = 1./M

    ## matrix definition
    mat = np.zeros((M-1,M-1))
    
    for i in range(M-2):
        mat[i,i] = 2.
        mat[i,i+1] = -1.
        mat[i+1,i] = -1.
    mat[M-2,M-2] = 2.
    
    mat = mat*(1./(h*h))
    #print("mat = ", mat)
    
    ## RHS
    X = np.linspace(0,1,M+1)
    RHS = np.zeros(M-1)
    for i in range(M-1):
        RHS[i] = rhs_fun(X[i+1])

    RHS[0] += A*(1./(h*h))
    RHS[M-2] += B*(1./(h*h))
    # print("RHS = ", RHS)

    ## solve the linear system
    sol = linalg.solve(mat,RHS)

    ## rebuild the solution vector
    Y = np.zeros(M+1)
    for i in range(1,M):
        Y[i] = sol[i-1]
    Y[0] = A
    Y[M] = B

    return Y

## numerical solution
Y = Dirichlet(M, f, alpha, beta)
X = np.linspace(0,1,M+1)
plt.plot(X,Y, marker="x", label = "solution obtenue numeriquement")


## exact solution
ex = np.zeros(M+1)
for i in range(0,M+1):
    ex[i] = sol_ex(X[i])

## plot
plt.plot(X,ex, label="solution exacte")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Probleme de Laplace avec conditions de Dirichlet")
plt.legend(loc="lower right")
plt.show()


#############  MIXED BC  ############
## matrix definition
mat = np.zeros((M-1,M-1))
h = 1/M

for i in range(M-2):
    mat[i,i] = 2.
    mat[i,i+1] = -1.
    mat[i+1,i] = -1.
mat[M-2,M-2] = 1.

mat = mat*(1./(h*h))

## RHS
RHS = np.zeros(M-1)
for i in range(M-1):
    RHS[i] = f_mixte(X[i+1])

## solve the linear system
sol = linalg.solve(mat,RHS)

## plot the solution
Y = np.zeros(M+1)
for i in range(1,M):
    Y[i] = sol[i-1]
Y[0] = 0.
Y[M] = Y[M-1]


plt.plot(X,Y, marker="x", label = "solution obtenue numeriquement")


## exact solution
ex = np.zeros(M+1)
for i in range(0,M+1):
    ex[i] = sol_mixte(X[i])

plt.plot(X,ex, label="solution exacte")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Probleme de Laplace avec conditions mixtes")
plt.legend(loc="lower right")
plt.show()


##############  CONVERGENCE RATE (Dirichlet) ###########

meshes = [25,50,100,200,400]
Err = []
liste_h = []
for M in meshes:

    ## mesh discretization
    X = np.linspace(0,1,M+1)

    ## numerical solution
    Y = Dirichlet(M, f, alpha, beta)

    ## exact solution
    ex = np.zeros(M+1)
    for i in range(0,M+1):
        ex[i] = sol_ex(X[i])

    ## error
    err = linalg.norm(Y - ex, np.inf)
    Err.append(err)
    liste_h.append(1/M)


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
