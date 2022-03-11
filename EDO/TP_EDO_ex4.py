import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#################  Exercice 4 : Etude du pendule
## etude du pendule
## dans l'ideal, il faudrait ecrire les resultats theoriques dans le sujet de TP

##### Parametres
g = 9.81
l = 0.1

def f(X,t):
    return np.array([X[1], -g/l*np.sin(X[0])])

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

## donnees simu
X0 = np.array([np.pi/3,0])
T = 2
N = 200

## remarque : pour T = 10 et N = 200, on voit tres bien l'effet d'amplification du schema d'Euler,
## on pourrait faire une question la-dessus

# Euler
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

## affichage de theta
plt.figure()
plt.plot(tps,theta_Euler,label="Euler")
plt.plot(tps,theta_milieu,label="point milieu")
plt.plot(tps,theta_ode,label="odeint")
plt.xlabel("t")
plt.ylabel("theta")
plt.title("Theta(t)")
plt.legend()
plt.show()

## affichage de omega
plt.figure()
plt.plot(tps,omega_Euler,label="Euler")
plt.plot(tps,omega_milieu,label="point milieu")
plt.plot(tps,omega_ode,label="odeint")
plt.xlabel("t")
plt.ylabel("omega")
plt.title("Omega(t)")
plt.legend()
plt.show()

## portrait de phases
plt.figure()
plt.plot(theta_Euler,omega_Euler,label="Euler")
plt.plot(theta_milieu,omega_milieu,label="point milieu")
plt.plot(theta_ode,omega_ode,label="odeint")
plt.xlabel("theta")
plt.ylabel("omega")
plt.title("Portrait de phases")
plt.legend()
plt.show()

## question importante : quel schema est utilise dans odeint ??

## ici ca serait sympa de pouvoir tracer un portrait de phases super detaille avec
## isoclines, fleches et isocontours donnes par 0.5*y^2 + a cos(x) = cste (cf cours)


## il faudrait rajouter un exo sur les schemas implicite, methode de Newton,
## et si possible mettre en evidence le fait que le schema n'est pas defini pour un pas de temps grand
