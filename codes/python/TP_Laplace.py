#################### IMPORTS ##############
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
##################  PARAMETERS  ###########

M = 100 # number of space sub-divisions


## RHS function
def f(x):
    #return x
    return x*x

## Exact solution
def sol_ex(x):
    #return -x**3/6. + x/6.
    return x/12. * (1.-x*x*x) + x/12.

def sol_mixte(x):
    return x/3. * (1.-x*x*x/4.)

##################  MAIN CODE  ############
## step
h = 1./M

## Dirichlet data
alpha = sol_ex(0)
beta = sol_ex(1)

## X = (x_j)
X = np.linspace(0,1,M+1)
# print("X = ", X)


#############  DIRICHLET BC  ############
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
RHS = np.zeros(M-1)
for i in range(M-1):
    RHS[i] = f(X[i+1])

RHS[0] += alpha*(1./(h*h))
RHS[M-2] += beta*(1./(h*h))
# print("RHS = ", RHS)

## solve the linear system
sol = linalg.solve(mat,RHS)

## plot the solution
Y = np.zeros(M+1)
for i in range(1,M):
    Y[i] = sol[i-1]
Y[0] = alpha
Y[M] = beta


#plt.plot(X,Y)


## exact solution
ex = np.zeros(M+1)
for i in range(0,M+1):
    ex[i] = sol_ex(X[i])

#plt.plot(X,ex)
#plt.show()


#############  MIXED BC  ############
## matrix definition
mat = np.zeros((M-1,M-1))

for i in range(M-2):
    mat[i,i] = 2.
    mat[i,i+1] = -1.
    mat[i+1,i] = -1.
mat[M-2,M-2] = 1.

mat = mat*(1./(h*h))
#print("mat = ", mat)

## RHS
RHS = np.zeros(M-1)
for i in range(M-1):
    RHS[i] = f(X[i+1])

#RHS[0] += alpha*(1./(h*h))
#RHS[M-2] += beta*(1./(h*h))
# print("RHS = ", RHS)

## solve the linear system
sol = linalg.solve(mat,RHS)

## plot the solution
Y = np.zeros(M+1)
for i in range(1,M):
    Y[i] = sol[i-1]
Y[0] = 0.
Y[M] = Y[M-1]


plt.plot(X,Y)


## exact solution
ex = np.zeros(M+1)
for i in range(0,M+1):
    ex[i] = sol_mixte(X[i])

plt.plot(X,ex)
plt.show()
