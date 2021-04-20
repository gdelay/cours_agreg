////////  Equation de transport avec conditions periodiques par differences finies sur (0,1)
clear

/////// Parametres
M = 160   // nb de subdivisions d'espace
// N = 120   // nb de pas de temps
a = 1.0  // vitesse
T = 0.5  // temps final


// alpha : donnee de Dirichlet
function z = alpha(t)
//    z = 1.0 + t
    z = 2.0
endfunction

// u0 : donnee initiale
function z = u0(x)
    z = 0.0 .* x + 1.0
endfunction

// sol_ex : solution exacte a l'instant final T
function z = sol_ex(x)
    [n,n1] = size(x)
    
    for i = 1:n
        if(x(i) < 0.5)
            z(i) = 2.0
        elseif(x(i) == 0.5)
                z(i) = 1.5
        else
            z(i) = 1.0
        end
    end
endfunction

//////////////// code principal

figure(1,"Figure_name",'equation de transport Dirichlet')


// calcul et affichage de la solution exacte
X = linspace(0,1,M+1)'
U_ex = sol_ex(X)
plot(X,U_ex,'o-b')

// boucle principale
CFL = [0.9,0.5,0.01]
col = ['*-black', 'x-red', '+-yellow']
for it = 1:3
c = CFL(it) // nombre de CFL

hx = 1.0/M // pas d'espace
ht = c*hx/a   // pas de temps
N = floor(T/ht)      // nombre de pas de temps

disp(N)

// sous-divisions du temps
TT = linspace(0,T,N+1)'


// initialisation des inconnues
U0 = u0(X)
U_am = U0
U_ce = U0

// boucle en temps
// schema decentre amont
for i = 1:N
    // stockage du pas de temps precedent
    U_temp = U_am

    // mise a jour des inconnues
    for j = 2:M+1
        U_am(j) = (1-c) * U_temp(j) + c * U_temp(j-1)
    end
    // condition de Dirichlet
    U_am(1) = alpha(i*ht)
end

// dernier pas de temps (de taille reduite)
c_bis=a*(T-N*ht)/hx
U_temp = U_am

// mise a jour des inconnues
for j = 2:M+1
    U_am(j) = (1-c_bis) * U_temp(j) + c_bis * U_temp(j-1)
end
// condition de Dirichlet
U_am(1) = alpha(T)


// affichage de la solution approchee par le schema decentre amont
//plot(X,U_am,'*-black')
plot(X,U_am,col(it))
end
// fin de la boucle principale

legend("solution exacte","c=0.9","c=0.5", "c=0.01")
xtitle("Solution pour T=0.5, M=160","Position (x)","Solution (u)")

//figure(2,"Figure_name",'etude de la convergence')
//plot(log10(vect_h),log10(vect_err),'o-b')
//xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
