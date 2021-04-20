////////  Equation de transport avec conditions periodiques par differences finies sur (0,1)
clear

/////// Parametres
M = 160   // nb de subdivisions d'espace
N = 120   // nb de pas de temps
a = 1.0  // vitesse
T = 0.5  // temps final

// u0 : donnee initiale
function z = u0(x)
    z = x.*x.*(1.0-x).*(1.0-x)
endfunction


//////////////// code principal


hx = 1.0/M // pas d'espace
ht = T/N   // pas de temps

// sous-divisions de l'espace et du temps
X = linspace(0,1,M+1)'
TT = linspace(0,T,N+1)'


// initialisation des inconnues
U0 = u0(X)
U_am = U0
U_ce = U0

// nombre de cfl
c = ht * a / hx

// boucle en temps
// schema decentre amont
for i = 1:N
    // stockage du pas de temps precedent
    U_temp = U_am

    // mise a jour des inconnues
    for j = 2:M+1
        U_am(j) = (1-c) * U_temp(j) + c * U_temp(j-1)
    end
    // condition periodique
    U_am(1) = U_am(M+1)
end
// schema centre
for i = 1:N
    // stockage du pas de temps precedent
    U_temp = U_ce

    // mise a jour des inconnues
    for j = 2:M
        U_ce(j) = -0.5*c*U_temp(j+1) + U_temp(j) + 0.5*c * U_temp(j-1)
    end
    U_ce(M+1) = -0.5*c*U_temp(2) + U_temp(M+1) + 0.5*c * U_temp(M)
    // condition periodique
    U_ce(1) = U_ce(M+1)
end

// donnee initiale
figure(1,"Figure_name",'equation de transport Dirichlet')
//plot(X,U0,'o-b')

// calcul et affichage de la solution exacte
pos_ex = X - a*T - floor(X-a*T)
U_ex = u0(pos_ex)
plot(X,U_ex,'o-b')


// affichage de la solution approchee par le schema decentre amont
plot(X,U_am,'*-black')

// affichage de la solution approchee par le schema centre
plot(X,U_ce,'x-r')

legend("solution exacte","schema decentre amont","schema centre")
xtitle("Solution pour T=0.5, N=120 et M=160","Position (x)","Solution (u)")

//figure(2,"Figure_name",'etude de la convergence')
//plot(log10(vect_h),log10(vect_err),'o-b')
//xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
