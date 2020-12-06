############################
##### IMPORTED MODULES #####
############################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import copy

################################
##### FUNCTION DEFINITIONS #####
################################
def make_grille(triangle):
    xmin=triangle['maillage']['xmin']
    xmax=triangle['maillage']['xmax']
    Nx=triangle['maillage']['Nx']
    return np.linspace(xmin,xmax,Nx)

def init(triangle):
    grille=triangle['maillage']['grille']
    w_max=triangle['init']['w_max']
    fct = 1-np.abs(grille)
    fct[fct<0]=0
    return fct

def upwind(triangle):
    w_0=triangle['calcul']['init']
    grille=triangle['maillage']['grille']
    nt=triangle['calcul']['nt']
    u_max=triangle['calcul']['u_max']
    w_max=triangle['init']['w_max']
    
    CFL=triangle['calcul']['CFL']
    adv=triangle['calcul']['adv']
    adv_lin = triangle['calcul']['adv_lin']
    
    delta_x=grille[1]-grille[0]
    
    w_n=w_0
    w_np1=np.zeros(w_0.shape)

    sol_list=list()
    sol_list.append(w_0)

    dt=CFL*delta_x/adv
    
    for t in range(nt):
        
        #+++++ Conditions limites +++++
        w_np1[0]=w_n[-2]
        w_np1[-1]=w_n[1]
        
        #+++++ Interieur du domaine +++++
        a_jp1=0.5*(w_n[2:]+w_n[1:-1])
        a_jm1=0.5*(w_n[1:-1]+w_n[0:-2])
        
        if adv_lin == 'True':
            w_np1[1:-1]=w_n[1:-1]-dt/delta_x/2*(adv+np.abs(adv))*(w_n[1:-1]-w_n[0:-2])
            
        elif adv_lin == 'False':
            w_np1[1:-1]=w_n[1:-1]-dt/delta_x/2*(a_jm1+np.abs(a_jm1))*(w_n[1:-1]-w_n[0:-2])-dt/delta_x/2*(a_jp1-np.abs(a_jp1))*(w_n[2:]-w_n[1:-1])
        else :
            print('Choix du type de calcul invalide')
            break
        
        #+++++ Mise a jour de la solution +++++
        w_n=copy.copy(w_np1)

        #+++++ Enregistrement de la solution +++++
        sol_list.append(w_n)

    return np.array(sol_list)


def visualisation(triangle):
    grille=triangle['maillage']['grille']
    sol_list=triangle['calcul']['resultat']
    w_0=triangle['calcul']['init']
    
    #----- Figures -----
    figure=plt.figure(figsize=[7.0,5.0])
    plt.xlim([triangle['maillage']['xmin'],triangle['maillage']['xmax']])
    #plt.ylim([0.2,1.4])
    plt.plot(grille,sol_list[-1],'ro', label ='Solution à t=50')
    plt.plot(grille,w_0,'k-', label='Solution à t=0')
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel(r'$w$')

    plt.savefig(triangle['visualisation']['Image']['Fichier'])
    plt.show()

#######################
##### SCRIPT PART #####
#######################
#***** Parametres *****
triangle=dict()

triangle['maillage']=dict()
triangle['init']=dict()
triangle['calcul']=dict()
triangle['visualisation']=dict()

#----- Maillage -----
triangle['maillage']['Nx']=101
triangle['maillage']['xmin']=-2.0
triangle['maillage']['xmax']=2.0
triangle['maillage']['grille']=make_grille(triangle)

#----- Solution initiale -----
triangle['init']['w_max']=1.0
triangle['calcul']['init']=init(triangle)

#----- Calcul -----
triangle['calcul']['adv']=1.0
triangle['calcul']['CFL']=0.5
triangle['calcul']['nt']=50
triangle['calcul']['u_max']=1.0
triangle['calcul']['adv_lin']='False'

#***** Schema numerique *****
triangle['calcul']['resultat']=upwind(triangle)

#***** Sortie graphique *****
triangle['visualisation']['Image']=dict()
triangle['visualisation']['Image']['Fichier']='Initiale_finale.png' # jpg, png, pdf, ps, eps et svg.

visualisation(triangle)
