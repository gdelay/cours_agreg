import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy.linalg as LA

#################  Exercice 4 : Schemas implicites
## on resout le systeme
## x' = x + 2y
## y' = 2x + y

############  Parametres
# conseil : mettre tous les parametres que vous pouvez etre amenes a modifier dans une section au debut du code
T = 1 # temps final
N = 10 # nbre de pas de temps

# fonction second membre f
def f(X,t):
    x = X[0] + 2*X[1]
    y = 2*X[0] + X[1]
    return np.array([x,y])
# gradient de f par rapport a l'espace
def df(X,t):
    return np.array([[1,2],[2,1]])
# condition initiale
X0 = np.array([5.,0.])

# solution exacte :
def sol(t):
    a = (X0[0]-X0[1])/2.
    b = (X0[0]+X0[1])/2.
    x = a * np.exp(-t) + b * np.exp(3*t)
    y = - a * np.exp(-t) + b * np.exp(3*t)
    return np.array([x,y])

###############  methodes de resolution numerique
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

## schema de Crank-Nicolson
## df : jacobienne de f par rapport a l'espace
def Crank_Nicolson(f,df,X0,N,T):
    dt = T/N # pas de temps
    t=0 # temps
    X = np.zeros([N+1,2])
    X[0] = X0
    ## boucle en temps
    for i in range(0,N):
        def F(z):
            return z - X[i] - 0.5 * dt * (f(X[i],t) + f(z,t+dt))
        def DF(z):
            return np.identity(2) - 0.5 * dt * df(z,t+dt)
        X[i+1] = Newton(F,DF,X[i],0.000001,10)
        t += dt
    return X


##################  affichage des resultats

sol_Euler = Euler_imp(f,df,X0, N, T)
sol_CN = Crank_Nicolson(f,df,X0, N, T)
x_Euler, y_Euler = sol_Euler[:,0] , sol_Euler[:,1]
x_CN, y_CN = sol_CN[:,0] , sol_CN[:,1]
t_Euler = np.linspace(0,1,N+1)

## solution exacte :
sol_exacte = sol(t_Euler)
sol_x , sol_y = sol_exacte[0,:] , sol_exacte[1,:]

## figure x(t) (Euler, CN et sol exacte)
plt.figure()
plt.plot(t_Euler,x_Euler,label="Euler implicite")
plt.plot(t_Euler,x_CN,label="CN")
plt.plot(t_Euler,sol_x,label="x(t)")
plt.xlabel("t")
plt.ylabel("x")
plt.legend()
plt.show()

## figure y(t) (Euler, CN et sol exacte)
plt.figure()
plt.plot(t_Euler,y_Euler,label="Euler implicite")
plt.plot(t_Euler,y_CN,label="CN")
plt.plot(t_Euler,sol_y,label="y(t)")
plt.xlabel("t")
plt.ylabel("y")
plt.legend()
plt.show()

## solution dans le plan des phases (Euler, CN et sol exacte)
plt.figure()
plt.plot(x_Euler,y_Euler,label="Euler : yn(xn) EI")
plt.plot(x_CN,y_CN,label="Euler : yn(xn) CN")
plt.plot(sol_x,sol_y,label="sol exacte : y(x)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()


################ calcul du taux de convergence
## calcul de l'erreur max commise par le schema d'Euler implicite
def erreur_Euler_imp(X0,N,T):
    sol_Euler = Euler_imp(f,df,X0, N, T)
    x_Euler , y_Euler = sol_Euler[:,0] , sol_Euler[:,1]
    tps = np.linspace(0,T,N+1)
    sol_exacte = np.transpose(sol(tps))
    sol_x , sol_y = sol_exacte[:,0] , sol_exacte[:,1]

    # erreur
    err_x = LA.norm(x_Euler - sol_x, np.inf)
    err_y = LA.norm(y_Euler - sol_y, np.inf)

    return max([err_x, err_y])

## calcul de l'erreur max commise par le schema de Crank-Nicolson
def erreur_CN(X0,N,T):
    sol_CN = Crank_Nicolson(f,df,X0, N, T)
    x_CN , y_CN = sol_CN[:,0] , sol_CN[:,1]
    tps = np.linspace(0,T,N+1)
    sol_exacte = np.transpose(sol(tps))
    sol_x , sol_y = sol_exacte[:,0] , sol_exacte[:,1]

    # erreur
    err_x = LA.norm(x_CN - sol_x, np.inf)
    err_y = LA.norm(y_CN - sol_y, np.inf)

    return max([err_x, err_y])


## erreur en fonction du pas de temps (echelle logarithmique)
pas = [50,100,200,500,1000]

err_Euler = []
err_CN = []
delta_t = []
for k in pas:
    err_Euler.append(erreur_Euler_imp(X0,k,T))
    err_CN.append(erreur_CN(X0,k,T))
    delta_t.append(T/k)

## plot de l'erreur en echelle logarithmique
log_delta_t = np.log10(delta_t)
log_err_Euler = np.log10(err_Euler)
log_err_CN = np.log10(err_CN)

plt.figure()
plt.plot(log_delta_t, log_err_Euler, label="Euler implicite")
plt.plot(log_delta_t, log_err_CN, label="Crank-Nicolson")
plt.xlabel("dt")
plt.ylabel("erreur")
plt.legend()
plt.show()

## droite proche de cette erreur
z_Euler = np.polyfit(log_delta_t,log_err_Euler,1)
print('z_Euler = ', z_Euler)
z_CN = np.polyfit(log_delta_t,log_err_CN,1)
print('z_Euler = ', z_CN)
