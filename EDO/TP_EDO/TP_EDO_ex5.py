import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#################  Exercice 5 : Portrait de phase et stabilite des points critiques
## systeme de Lotka - Volterra


T = 30
N = 1000

## fonction second membre
def f(X,t):
    return np.array([X[0]*(2.-X[1]/10), X[1]*(X[0]/10 - 4.)])

X0 = np.array([100., 20.])
t = np.linspace(0,T,N+1)

## resolution
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
plt.plot(pt_critique_x, pt_critique_y, marker = '*', label="points critiques",linestyle='None')
# isoclines
plt.plot(isocline1_x, isocline1_y, label="isocline1")
plt.plot(isocline2_x, isocline2_y, label="isocline2")
## autres orbites
for k in np.linspace(0.1,2.0,10):
    X0 = np.array([100.*k, 20.*k])
    sol = integrate.odeint(f,X0,t)
    x, y = sol[:,0], sol[:,1]
    plt.plot(x,y)
## ajout de f
x=np.linspace(3,200,30)
y=np.linspace(3,200,30)
X1,Y1=np.meshgrid(x,y)
DX1,DY1=f([X1,Y1],0.0)
M=np.hypot(DX1,DY1)
M[M==0]=1.
DX1/= M
DY1/= M
plt.quiver(X1,Y1,DX1,DY1,M)
## finalisation
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()
