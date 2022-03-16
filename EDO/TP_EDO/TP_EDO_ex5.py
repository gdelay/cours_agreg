import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#################  Exercice 3 : Portrait de phase et stabilite des points critiques
## systeme de Lotka - Volterra
## dans l'ideal, il faudrait ecrire les resultats theoriques dans le sujet de TP


T = 30
N = 1000

## fonction second membre
def f(X,t):
    return np.array([X[0]*(2.-X[1]/10), X[1]*(X[0]/10 - 4.)])

X0 = np.array([100., 20.])
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
pt_critique_x = np.array([0. , 40.])
pt_critique_y = np.array([0. , 20.])

## isoclines
isocline1_x = np.array([0.,200.])
isocline1_y = np.array([20.,20.])

isocline2_x = np.array([40.,40.])
isocline2_y = np.array([0.,200.])

## affichage de la solution dans le plan des phases
plt.figure()
plt.plot(x,y,label="solution")
plt.plot(pt_critique_x, pt_critique_y, marker = '*', label="point critique")
plt.plot(isocline1_x, isocline1_y, marker = '*', label="isocline1")
plt.plot(isocline2_x, isocline2_y, marker = '*', label="isocline2")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()


## on peut ameliorer cette partie en utilisant ce qui est propose sur cette page
## (pour tracer le champ de f)
## https://les-mathematiques.net/vanilla/index.php?p=discussion/1255755#Comment_1255755
