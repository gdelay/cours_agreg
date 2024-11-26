import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy.linalg as LA

#################  Exercice 2 : Etude du pendule, simulations numeriques
## etude du pendule
##### Parametres
g = 9.81
l = 5.
X0 = np.array([np.pi/3,0])
T = 12
N = 200

def f(X,t):
    return np.array([X[1], -g/l*np.sin(X[0])])
def df(X,t):
    return np.array([[0, 1],[-g/l*np.cos(X[0]), 0]])

######## schemas
## schema d'Euler explicite
def Euler(X0,N,T):
    dt = T/N
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        print("iter = ", i)
        F = np.array([X[i,1],-g/l*np.sin(X[i,0])])
        X[i+1] = X[i] + dt * F

    return X

## schema du point milieu
## est-ce que ca ne serait pas mieux de mettre ce schema dans un autre exo plutot qu'ici ??
def milieu(X0,N,T):
    dt = T/N
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        print("iter = ", i)
        F = np.array([X[i,1],-g/l*np.sin(X[i,0])])
        X2 = X[i] + 0.5 * dt * F
        F2 = np.array([X2[1],-g/l*np.sin(X2[0])])
        X[i+1] = X[i] + dt * F2

    return X

## methode de Newton
def Newton(f,df,x0,tol,maxiter):
    x = x0
    it = 0
    delta_norm = 10 * tol
    while delta_norm > tol and it < maxiter :
        delta_x = np.dot( LA.inv( df(x) ) , f(x) )
        x = x - delta_x
        it += 1
        delta_norm = LA.norm(delta_x)
    return x

## schema d'Euler implicite
## df : jacobienne de f par rapport a l'espace
def Euler_imp(f,df,X0,N,T):
    dt = T/N # pas de temps
    t=0 # temps
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        def F(z):
            return z - X[i] - dt * f(z,t+dt)
        def DF(z):
            return np.identity(2) - dt * df(z,t+dt)
        X[i+1] = Newton(F,DF,X[i],0.000001,10)
        t += dt
    return X


# Euler explicite
X_Euler = Euler(X0,N,T)
theta_Euler , omega_Euler = X_Euler[:,0] , X_Euler[:,1]
# point milieu
X_milieu = milieu(X0,N,T)
theta_milieu , omega_milieu = X_milieu[:,0] , X_milieu[:,1]
# temps
tps = np.linspace(0,T,N+1)
# odeint
sol_ODE = integrate.odeint(f,X0,tps)
theta_ode , omega_ode = sol_ODE[:,0] , sol_ODE[:,1]
# Euler implicite
sol_EI = Euler_imp(f,df,X0,N,T)
theta_EI, omega_EI = sol_EI[:,0], sol_EI[:,1]


## affichage de theta
plt.figure()
plt.plot(tps,theta_Euler,label="Euler explicite")
plt.plot(tps,theta_milieu,label="point milieu")
plt.plot(tps,theta_ode,label="odeint")
plt.plot(tps,theta_EI,label="Euler implicite")
plt.xlabel("t")
plt.ylabel("theta")
plt.title("Theta(t)")
plt.legend()
plt.show()

## affichage de omega
plt.figure()
plt.plot(tps,omega_Euler,label="Euler explicite")
plt.plot(tps,omega_milieu,label="point milieu")
plt.plot(tps,omega_ode,label="odeint")
plt.plot(tps,omega_EI,label="Euler implicite")
plt.xlabel("t")
plt.ylabel("omega")
plt.title("Omega(t)")
plt.legend()
plt.show()

## portrait de phases
plt.figure()
plt.plot(theta_Euler,omega_Euler,label="Euler explicite")
plt.plot(theta_milieu,omega_milieu,label="point milieu")
plt.plot(theta_ode,omega_ode,label="odeint")
plt.plot(theta_EI,omega_EI,label="Euler implicite")
plt.xlabel("theta")
plt.ylabel("omega")
plt.title("Comportement des differents schemas dans l'espace des phases")
plt.legend()
plt.show()

