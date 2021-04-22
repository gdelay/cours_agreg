#################### IMPORTS ##############
import numpy as np
import matplotlib.pyplot as plt
##################  PARAMETERS  ###########

N = 100 # number of time  sub-divisions
M = 100 # number of space sub-divisions

T = 1.  # final time
a = 1.  # transport velocity

# initial solution
def u0(x):
    # return x
    return x**2


##################  MAIN CODE  ############

## steps
ht = T/N
hx = 1./M

## X = (x_j)
X = np.linspace(0,1,M+1)

################  SCHEMA CENTRE  ##########

## initialisation of the solution
U = np.zeros(M+1)
for j in range(M+1):
    U[j] = u0(X[j])

## time loop
CFL = a * ht / (2*hx)

for n in range(N):
    oldj = U[0]
    for j in range(1,M):
        oldjj = oldj
        oldj  = U[j]
        U[j] -= CFL * (U[j+1] - oldjj)
    U[M] -= CFL * (U[M] - oldj)
## plot final solution
plt.plot(X,U)
plt.show()


################  SCHEMA DECENTRE  ##########

## initialisation of the solution
U = np.zeros(M+1)
for j in range(M+1):
    U[j] = u0(X[j])

## time loop
CFL = a * ht / hx
print("CFL = ", CFL)

for n in range(N):
    oldj = U[0]
    for j in range(1,M+1):
        oldjj = oldj
        oldj  = U[j]
        U[j] -= CFL * (U[j] - oldjj)

## plot final solution
plt.plot(X,U)
plt.show()



