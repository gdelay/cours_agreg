import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#################  Exercice 2 : Portrait de phase et stabilite des points critiques
## peut etre que ce serait bien de remplacer cet exercice par quelque chose qui a deja ete fait en cours comme par exemple le portrait de phase du pendule (cf exo 4)
##

T = 10
N = 1000

## fonction second membre
def f(X,t):
    return np.array([-2*X[0]-X[1]+6, 3*X[0]-4*X[1]-7])

X0 = np.array([3., 1.])
t = np.linspace(0,T,N+1)

## resolution auto (quel schema ??)
sol = integrate.odeint(f,X0,t)

## affichage du resultat en fonction du temps
x , y = sol[:,0], sol[:,1]
plt.figure()
plt.plot(t,x,label="x(t)")
plt.plot(t,y,label="y(t)")
plt.xlabel("t")
plt.ylabel("x ou y")
plt.legend()
plt.show()

## point critique du systeme
pt_critique = np.array([31./11. , 4./11.])

## affichage de la solution dans le plan des phases
plt.figure()
plt.plot(x,y,label="solution")
plt.plot(pt_critique[0], pt_critique[1], marker = '*', label="point critique")
#plt.scatter(pt_critique[0], pt_critique[1], marker = '*',label="point critique")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()
