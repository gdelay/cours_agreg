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
N = 10 # nbre de pas de temps

# fonction second membre f
def f(t,X):
    x = X[0] + 2*X[1]
    y = 2*X[0] + X[1]
    return np.array([x,y])
# condition initiale
X0 = np.array([0.,1.])

############ Resolution
dt = T/N # pas de temps
print("dt = ", dt)

# solution X
X = np.zeros([N+1,2])
X[0] = X0
t = 0

## boucle en temps
for i in range(0,N):
    print("iter = ",i)
    F = f(t,X[i])
    X[i+1] = X[i] + dt * F
    t += dt

x_Euler , y_Euler = X[:,0] , X[:,1]

print(X)

## affichage du resultat
plt.figure()
plt.plot(x_Euler,y_Euler,label="solution")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()


## TODO : have a look to this website :
#  https://numerical-analysis.readthedocs.io/en/latest/ODE/ODE.html
