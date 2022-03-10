import numpy as np
import matplotlib.pyplot as plt

######### Exercice 1 : schema d'Euler explicite
## on resout le systeme
## x' = x + 2y
## y' = 2x + y
##

############  Parametres
# conseil : mettre tous les parametres que vous pouvez etre amenes a modifier dans une section au debut du code
T = 1.0 # temps final
N = 100 # nbre de pas de temps

# fonction second membre f
def f(t,X):
    x = X[0] + 2*X[1]
    y = 2*X[0] + X[1]
    return np.array([x,y])
# condition initiale
X0 = np.array([5.,0.])

# solution exacte :
def sol(t):
    a = (X0[0]-X0[1])/2.
    b = (X0[0]+X0[1])/2.
    x = a * np.exp(-t) + b * np.exp(3*t)
    y = - a * np.exp(-t) + b * np.exp(3*t)
    return np.array([x,y])

############ Resolution
dt = T/N # pas de temps
print("dt = ", dt)


## Euler explicite
def Euler_exp(X0,N,dt):
    t=0
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        print("iter = ", i)
        F = f(t,X[i])
        X[i+1] = X[i] + dt * F
        t += dt
    return X

sol_Euler = Euler_exp(X0, N, dt)
x_Euler, y_Euler = sol_Euler[:,0] , sol_Euler[:,1]
t_Euler = np.linspace(0,1,N+1)

## affichage du resultat (X,Y)
plt.figure()
plt.plot(x_Euler,y_Euler,label="solution y(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()

## affichage du resultat (t,X)
plt.figure()
plt.plot(t_Euler,x_Euler,label="x(t)")
plt.xlabel("t")
plt.ylabel("x")
plt.legend()
plt.show()

## affichage du resultat (t,Y)
plt.figure()
plt.plot(t_Euler,y_Euler,label="y(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.show()

## solution exacte :
sol_x = []
sol_y = np.zeros(N+1)
for i in range(0,N+1):
    t = i*dt
    sol_xi, sol_yi = sol(t)
    sol_x.append(sol_xi)
    sol_y[i] = sol_yi


plt.figure()
plt.plot(t_Euler,y_Euler,label="yn")
plt.plot(t_Euler,sol_y,label="y(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.show()


## faire etude du taux de convergence




## ce serait bien aussi de reorganiser cette partie





## TODO : have a look to this website :
#  https://numerical-analysis.readthedocs.io/en/latest/ODE/ODE.html
