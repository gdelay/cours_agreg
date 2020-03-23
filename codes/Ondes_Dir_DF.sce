//////// Equation de la chaleur avec conditions de Dirichlet (par differences finies sur (0,1))
clear

/////// Parametres
M = 50   // nb de subdivisions en espace
N = 400   // nb de subdivisions en temps
c = 0.2     // vitesse onde
T = 0.4   // temps final

N_int = 200 // sortie intermediaire

// u0 : donnee initiale
function z = u0(x)
    for i = 1:length(x)
       if(x(i)<=0.3 || x(i) >= 0.7)
            z(i) = 0
        elseif(x(i) < 0.5)
            z(i) = x(i)-0.3
        else
            z(i) = 0.7-x(i)
        end
    end
endfunction

// v0 : vitesse initiale
function z = v0(x)
    for i=1:length(x)
        z(i) = 0.0
    end
endfunction

// f : terme source
function z = f(x)
    z=0
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
U0 = u0(X)   // solution initiale
V0 = v0(X)   // vitesse initiale
U_prec(1:M-1) = U0(2:M)
U(1:M-1) = U0(2:M) + ht .* V0(2:M)


// nombre de cfl
CFL = c*ht/hx

// boucle en temps
for i = 2:N
    // stockage du pas de temps precedent
    U_temp = U

    // mise a jour de la solution
    U(1) = 2*U_temp(1) - U_prec(1) + CFL*CFL*(U_temp(2) - 2 * U_temp(1))// + ht*ht * F(1)
    for j = 2:M-2
        U(j) = 2*U_temp(j) - U_prec(j) + CFL*CFL * (U_temp(j+1) -2* U_temp(j) + U_temp(j-1)) //+ ht*ht * F(j)
    end
    U(M-1) = 2*U_temp(M-1)-U_prec(M-1) + CFL*CFL * (-2*U_temp(M-1)+U_temp(M-2)) //+ ht*ht * F(M-1)

    // mise a jour de U_prec
    U_prec = U_temp

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
figure(1,"Figure_name",'conditions de Dirichlet homogenes')
plot(X,U0,'o-b')

// affichage de la solution intermediaire
plot(X,U_int,'x-r')

// affichage de la solution finale
plot(X,U_fin,'*-black')

xtitle("Equation des ondes","position (x)","solution (u)")
legend("t=0","t=0.2","t=0.4")

//figure(2,"Figure_name",'etude de la convergence')
//plot(log10(vect_h),log10(vect_err),'o-b')
//xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
