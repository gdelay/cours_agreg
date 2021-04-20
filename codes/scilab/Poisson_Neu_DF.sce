////////  Probleme de Poisson-Neumann par differences finies sur (0,1)
clear

/////// Parametres
M = 100   // nb de subdivisions
h = 1.0/M // pas
k = 1.0    // raideur


vect_M = [10,20,50,100,200]
//vect_M=[M]
for l=1:length(vect_M)
M = vect_M(l)
h=1.0/M

// sous-divisions de l'espace
X = linspace(0,1,M+1)

// f : terme source
function z = f(x)
    z=x-0.5
endfunction

//function z = f(x)
//    if(x < 0.5)
//        z = 1.0
//    elseif(x == 0.5)
//        z = 0.0
//    else
//        z = -1.0
//    end
//endfunction

// u : solution exacte
function z = u(x)
    z = 0.25*x.*x - (1.0/6.0) * x.*x.*x - 1.0/24.0
endfunction

//function z = u(x)
//    for i=1:length(x)
//        if(x(i) < 0.5)
//            z(i) = -0.5*x(i)*x(i) + 1.0/8.0
//        else
//            z(i) = 0.5*x(i)*x(i) -x(i) +3.0/8.0
//        end
//    end
//endfunction

// sous-divisions de l'espace
X = linspace(0,1,M+1)'


// calcul de la matrice de Poisson Dirichlet
A = zeros(M,M)
A(1,1) = 1
A(1,2) = -1
A(1,M) = h*h
for i=2:M-2
    A(i,i-1) = -1
    A(i,i) = 2
    A(i,i+1) = -1
    A(i,M) = h*h
end
A(M-1,M-2) = -1
A(M-1,M-1) = 1
A(M-1,M) = h*h

for i=1:M-1
    A(M,i) = h*h
end

A = k.*A./(h*h)

// calcul du second membre F
for i=1:M-1
    F(i) = f(X(i+1))
end
F(M) = 0

// calcul de la solution num U
U_temp = A \ F
U(1) = U_temp(1)
U(2:M) = U_temp(1:M-1)
U(M+1) = U_temp(M-1)
alpha = U_temp(M)

// affichage de la solution exacte
Uex = u(X)
figure(1,"Figure_name",'probleme de Dirichlet homogene')
plot(X,Uex,'o-b')

// affichage de la solution approchee
plot(X,U,'o-r')

// calcul de l'erreur
err = abs(Uex - U)
disp(max(err))
disp(alpha)
vect_err(l) = max(err)
vect_h(l) = h
end

figure(2,"Figure_name",'etude de la convergence')
plot(log10(vect_h),log10(vect_err),'o-b')
xtitle("Erreur du schema en fonction de h (echelle logarithmique)","Pas (h)","Erreur")
//legend("k=5","k=1","k=0.5")
//xtitle("Solution du problème de Poisson-Dirichlet","Position (x)","Solution (u)")
