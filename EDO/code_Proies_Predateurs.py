###################################################################
# Système Proies-Prédateurs (Modèle de Lokta-Volterra)
###################################################################

###################################################################
# Importations
###################################################################

import numpy as np                
 
from scipy.integrate import odeint
#from scipy.optimize import fsolve

import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches

###################################################################
# Fonctions
###################################################################

def Lokta_Volterra():
    def F(v,t):   # Prend v=[x,y] et renvoie [a x y - b x, -c x y + d y]   
        return np.array([Coeff[0]*v[0]*v[1]-Coeff[1]*v[0], -Coeff[2]*v[0]*v[1]+Coeff[3]*v[1]])
    
    Coeff=[1.0, 2.0, 1.0, 3.0]
    tmin, tmax, tpas = 0.0, 10.0, 0.01  
    t = np.arange(tmin,tmax,tpas)    
    
    Sol0=np.array([1.0, 1.0])
    Sol=odeint(F,Sol0,t) 
    
    plt.figure()
    plt.plot(t,Sol[:,0],'r',t,Sol[:,1],'b',linewidth=3)
    plt.title('Evolution des populations',size=40,color='m',fontweight='bold')
    label1 = mpatches.Patch(color='b', label='Proies')
    label2 = mpatches.Patch(color='r', label='Prédateurs')
    lg=plt.legend(loc='upper right', fontsize='xx-large', shadow=True, fancybox = True, handles=[label1, label2])
    lg.get_frame().set_facecolor('cyan')
    lg.get_frame().set_edgecolor('orchid')
    lg.get_frame().set_linewidth(5.0)  
    
    plt.figure()
    plt.plot(Sol[:,0],Sol[:,1],linewidth=3)
    plt.title('Diagramme de phase',size=40,color='m',fontweight='bold')
    plt.xlabel('Nombre de prédateurs',size=40,color='r',fontweight='bold')
    plt.ylabel('Nombre de proies',size=40,color='b',fontweight='bold')
    
    plt.figure()
    plt.grid()
    for k in np.linspace(0.2,1.0,10):
        v=odeint(F,k*Sol0,t)
        plt.plot(v[:,0],v[:,1],linewidth=3)
    x=np.linspace(0,18,20)
    y=np.linspace(0,16,20)
    plt.title('Diagramme de phase et champ de vecteurs',size=40,color='m',fontweight='bold')
    plt.xlabel('Nombre de prédateurs',size=40,color='r',fontweight='bold')
    plt.ylabel('Nombre de proies',size=40,color='b',fontweight='bold')
    X1,Y1=np.meshgrid(x,y)
    DX1,DY1=F([X1,Y1],0.0)
    M=np.hypot(DX1,DY1)
    M[M==0]=1.
    DX1/= M
    DY1/= M
    plt.quiver(X1,Y1,DX1,DY1,M)
    plt.show()
    
###################################################################
# Script principal de test
###################################################################

plt.close('all')

Lokta_Volterra()

###################################################################







