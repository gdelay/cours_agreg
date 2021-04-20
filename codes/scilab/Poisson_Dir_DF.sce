////////  Probleme de Poisson Dirichlet par differences finies sur (0,1)
clear

/////// Parametres
M = 100   // nb de subdivisions
h = 1.0/M // pas
k = 1.0    // raideur


vect_M = [10,20,50,100,200]
for l=1:length(vect_M)
M = vect_M(l)
h=1.0/M

// sous-divisions de l'espace
X = linspace(0,1,M+1)

// f : terme source
function z = f(x)
    //z = 1.0
    //z = x
    z=x.*x
endfunction

// u : solution exacte
function z = u(x)
    // z = x.*(1.0-x) ./ (2.0*k)
    //z = x.*(1.0-x.*x) ./ (6.0*k)
    z = x.*(1.0-x.*x.*x) / (12.0*k)
endfunction


// sous-divisions de l'espace
X = linspace(0,1,M+1)'


// calcul de la matrice de Poisson Dirichlet
A = zeros(M-1,M-1)
A(1,1) = 2
A(1,2) = -1
for i=2:M-2
    A(i,i-1) = -1
    A(i,i) = 2
    A(i,i+1) = -1
end
A(M-1,M-2) = -1
A(M-1,M-1) = 2

A = k.*A./(h*h)

// calcul du second membre F
for i=1:M-1
    F(i) = f(X(i+1))
end

// calcul de la solution num U
U_temp = A \ F
U(1) = 0
U(2:M) = U_temp
U(M+1) = 0

// affichage de la solution exacte
Uex = u(X)
figure(1,"Figure_name",'probleme de Dirichlet homogene')
plot(X,Uex,'o-b')

// affichage de la solution approchee
plot(X,U,'o-r')

// calcul de l'erreur
err = abs(Uex - U)
disp(max(err))
vect_err(l) = max(err)
vect_h(l) = h
end

figure(2,"Figure_name",'etude de la convergence')
plot(log10(vect_h),log10(vect_err),'o-b')
xtitle("Erreur du schema en fonction de h (echelle logarithmique)","log(h)","log(erreur)")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
