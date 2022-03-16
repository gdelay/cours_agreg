import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

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

############ Resolution par le schema d'Euler explicite (et pt milieu)

## Euler explicite
def Euler_exp(X0,N,T):
    dt = T/N # pas de temps
    t=0 # temps
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        F = f(t,X[i])
        X[i+1] = X[i] + dt * F
        t += dt
    return X

## point milieu
def pt_milieu(X0,N,T):
    dt = T/N # pas de temps
    t=0 # temps
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        F1 = f(t,X[i])
        X1 = X[i] + 0.5*dt * F1
        F = f(t+0.5*dt, X1)
        X[i+1] = X[i] + dt * F
        t += dt
    return X

##################  affichage des resultats

sol_Euler = Euler_exp(X0, N, T)
x_Euler, y_Euler = sol_Euler[:,0] , sol_Euler[:,1]
t_Euler = np.linspace(0,1,N+1)

## solution exacte :
sol_exacte = sol(t_Euler)
sol_x , sol_y = sol_exacte[0,:] , sol_exacte[1,:]

## figure x(t) (Euler et sol exacte)
plt.figure()
plt.plot(t_Euler,x_Euler,label="xn")
plt.plot(t_Euler,sol_x,label="x(t)")
plt.xlabel("t")
plt.ylabel("x")
plt.legend()
plt.show()

## figure y(t) (Euler et sol exacte)
plt.figure()
plt.plot(t_Euler,y_Euler,label="yn")
plt.plot(t_Euler,sol_y,label="y(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.show()

## solution dans le plan des phases (Euler et sol exacte)
plt.figure()
plt.plot(x_Euler,y_Euler,label="Euler : yn(xn)")
plt.plot(sol_x,sol_y,label="sol exacte : y(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()

################ calcul du taux de convergence
## calcul de l'erreur max commise par le schema d'Euler
def erreur_Euler_exp(X0,N,T):
    sol_Euler = Euler_exp(X0, N, T)
    x_Euler , y_Euler = sol_Euler[:,0] , sol_Euler[:,1]
    tps = np.linspace(0,T,N+1)
    sol_exacte = np.transpose(sol(tps))
    sol_x , sol_y = sol_exacte[:,0] , sol_exacte[:,1]

    # erreur
    err_x = LA.norm(x_Euler - sol_x, np.inf)
    err_y = LA.norm(y_Euler - sol_y, np.inf)

    return max([err_x, err_y])

## calcul de l'erreur max commise par le schema du point milieu
def erreur_milieu(X0,N,T):
    sol_mil = pt_milieu(X0, N, T)
    x_mil , y_mil = sol_mil[:,0] , sol_mil[:,1]
    tps = np.linspace(0,T,N+1)
    sol_exacte = np.transpose(sol(tps))
    sol_x , sol_y = sol_exacte[:,0] , sol_exacte[:,1]

    # erreur
    err_x = LA.norm(x_mil - sol_x, np.inf)
    err_y = LA.norm(y_mil - sol_y, np.inf)

    return max([err_x, err_y])


## erreur en fonction du pas de temps (echelle logarithmique)
pas = [50,100,200,500,1000]

err_Euler = []
err_milieu = []
delta_t = []
for k in pas:
    err_Euler.append(erreur_Euler_exp(X0,k,T))
    err_milieu.append(erreur_milieu(X0,k,T))
    delta_t.append(T/k)

## plot de l'erreur en echelle logarithmique
log_delta_t = np.log10(delta_t)
log_err_Euler = np.log10(err_Euler)
log_err_milieu = np.log10(err_milieu)

plt.figure()
plt.plot(log_delta_t, log_err_Euler, label="Euler")
plt.plot(log_delta_t, log_err_milieu, label="Point milieu")
plt.xlabel("dt")
plt.ylabel("erreur")
plt.legend()
plt.show()

## droite proche de cette erreur
z_Euler = np.polyfit(log_delta_t,log_err_Euler,1)
print('z_Euler = ', z_Euler)
z_milieu = np.polyfit(log_delta_t,log_err_milieu,1)
print('z_Euler = ', z_milieu)
