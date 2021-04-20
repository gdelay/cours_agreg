////////  Equation de transport avec conditions de Dirichlet par differences finies sur (0,1)
clear

/////// Parametres
M = 20   // nb de subdivisions d'espace
N = 15   // nb de pas de temps
a = 1.5  // vitesse
T = 0.5  // temps final

// u0 : donnee initiale
function z = u0(x)
    z = 1.0-x
endfunction


//////////////// code principal


hx = 1.0/M // pas d'espace
ht = T/N   // pas de temps

// sous-divisions de l'espace et du temps
X = linspace(0,1,M+1)'
TT = linspace(0,T,N+1)'


// initialisation des inconnues
U = u0(X)
U0 = U

// nombre de cfl
c = ht * a / hx

// boucle en temps
for i = 1:N
    // stockage du pas de temps precedent
    U_temp = U

    // mise a jour des inconnues
    for j = 2:M+1
        U(j) = (1-c) * U_temp(j) + c * U_temp(j-1)
    end
end

// donnee initiale
figure(1,"Figure_name",'equation de transport Dirichlet')
//plot(X,U0,'o-b')

// affichage de la solution approchee
plot(X,U,'*-black')



legend("a=0","a=0.5","a=1","a=1.5")

//figure(2,"Figure_name",'etude de la convergence')
//plot(log10(vect_h),log10(vect_err),'o-b')
//xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
