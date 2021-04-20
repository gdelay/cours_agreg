//////// Equation de la chaleur avec conditions de Dirichlet (par differences finies sur (0,1))
clear

/////// Parametres
M = 100   // nb de subdivisions en espace
N = 4000   // nb de subdivisions en temps
nu = 0.5  // diffusion
T = 0.1   // temps final

N_int = 2000 // sortie intermediaire

// u0 : donnee initiale
function z = u0(x)
    z=x.*(1-x)
endfunction

// f : terme source
function z = f(x)
    z=-1
endfunction

///////// code principal

hx = 1.0/M // pas d'espace
ht = T/N   // pas de temps

// sous-divisions de l'espace
X = linspace(0,1,M+1)'


// calcul du second membre F
for i=1:M-1
    F(i) = f(X(i+1))
end

// initialisation de la solution
U0 = u0(X)
U(1:M-1) = U0(2:M)

// nombre de cfl
c = nu*ht/(hx*hx)

// boucle en temps
for i = 1:N
    // stockage du pas de temps precedent
    U_temp = U

    // mise a jour de la solution
    U(1) = (1-2*c)*U_temp(1) + c * U_temp(2) + ht * F(1)
    for j = 2:M-2
        U(j) = c * U_temp(j-1) + (1-2*c) * U_temp(j) + c*U_temp(j+1) + ht * F(j)
    end
    U(M-1) = c*U_temp(M-2) + (1-2*c) * U_temp(M-1) + ht * F(M-1)

    // sortie intermediaire
    if(i == N_int)
        U_in = U
    end
end
U_fin(1) = 0.0
U_fin(2:M) = U(1:M-1)
U_fin(M+1) = 0.0

U_int(1) = 0.0
U_int(2:M) = U_in
U_int(M+1) = 0.0

// donnee initiale
figure(1,"Figure_name",'probleme de Dirichlet homogene')
plot(X,U0,'o-b')

// affichage de la solution intermediaire
plot(X,U_int,'x-r')

// affichage de la solution finale
plot(X,U_fin,'*-black')

xtitle("Equation de la chaleur","position (x)","solution (u)")
legend("t=0","t=0.05","t=0.1")

//figure(2,"Figure_name",'etude de la convergence')
//plot(log10(vect_h),log10(vect_err),'o-b')
//xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
