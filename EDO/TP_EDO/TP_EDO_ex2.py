import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#################  Exercice 2 : Portrait de phase et stabilite des points critiques
## peut etre que ce serait bien de remplacer cet exercice par quelque chose qui a deja ete fait en cours comme par exemple le portrait de phase du pendule (cf exo 4)
##
##### Parametres
g = 9.81
l = 5
T = 12
N = 1000

## fonction second membre
def f(X,t):
    return np.array([X[1], -g/l*np.sin(X[0])])


tps = np.linspace(0,T,N+1)



###### portrait de phase
plt.figure()
tps = np.linspace(0,T,N+1)
## premier noeud (0)
for k in np.linspace(0.2,3.,10):
    X0 = np.array([k,0])
    sol_ODE = integrate.odeint(f,X0,tps)
    theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
    plt.plot(theta_ode,omega_ode)
## deuxieme noeud (2pi)
for k in np.linspace(0.2,3.,10):
    X0 = np.array([2*np.pi+k,0])
    sol_ODE = integrate.odeint(f,X0,tps)
    theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
    plt.plot(theta_ode,omega_ode)
## niveaux superieurs d'energie
## au-dessus
T = 4
tps = np.linspace(0,T,N+1)
for k in np.linspace(0.2,3.,10):
    X0 = np.array([-np.pi,k])
    sol_ODE = integrate.odeint(f,X0,tps)
    theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
    plt.plot(theta_ode,omega_ode)
## au-dessous
for k in np.linspace(0.2,3.,10):
    X0 = np.array([3*np.pi,-k])
    sol_ODE = integrate.odeint(f,X0,tps)
    theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
    plt.plot(theta_ode,omega_ode)
## points critiques
pt_critique_x = np.array([-np.pi, 0., np.pi, 2*np.pi, 3*np.pi])
pt_critique_y = np.array([0., 0., 0.,0., 0.])


plt.plot(pt_critique_x, pt_critique_y, marker = '*', label="points critiques",linestyle='None')
plt.xlabel("theta")
plt.ylabel("omega")
plt.title("Portrait de phases")
plt.legend()
plt.show()



#################### zoom sur (0,0)
T = 12
tps = np.linspace(0,T,N+1)
plt.figure()
## premier noeud (0)
for k in np.linspace(0.2,3.,10):
    X0 = np.array([k,0])
    sol_ODE = integrate.odeint(f,X0,tps)
    theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
    plt.plot(theta_ode,omega_ode)
## champ de f
x=np.linspace(-3.,3.,15)
y=np.linspace(-3,3,15)
X1,Y1=np.meshgrid(x,y)
DX1,DY1=f([X1,Y1],0.0)
M=np.hypot(DX1,DY1)
M[M==0]=1.
DX1/= M
DY1/= M
plt.quiver(X1,Y1,DX1,DY1,M)
## isoclines
x0 = np.array([0,0])
y0 = np.array([-np.pi,np.pi])
plt.plot(x0,y0,label="isocline en y")
x0 = np.array([-np.pi,np.pi])
y0 = np.array([0,0])
plt.plot(x0,y0,label="isocline en x")

#points critiques
pt_critique_x = np.array([-np.pi, 0., np.pi])
pt_critique_y = np.array([0., 0., 0.])
plt.plot(pt_critique_x, pt_critique_y, marker = '*', label="points critiques",linestyle='None')

plt.xlabel("theta")
plt.ylabel("omega")
plt.title("Portrait de phases (zoom sur (0,0))")
plt.axis('scaled')
plt.legend()
plt.show()
