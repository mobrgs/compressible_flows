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
def make_grille(step):
    xmin=step['maillage']['xmin']
    xmax=step['maillage']['xmax']
    Nx=step['maillage']['Nx']
    return np.linspace(xmin,xmax,Nx)

def init(step):
    grille=step['maillage']['grille']
    rho_max=step['init']['rho_max']
    fct = np.zeros(np.size(grille))
    fct[grille<0]=rho_max
    fct[grille>=0]=rho_max/2
    return fct

def upwind(step):
    rho_0=step['calcul']['init']
    grille=step['maillage']['grille']
    nt=step['calcul']['nt']
    u_max=step['calcul']['u_max']
    rho_max=step['init']['rho_max']
    
    CFL=step['calcul']['CFL']
    adv=step['calcul']['adv']
    stationnaire = step['calcul']['stationnaire']
    
    delta_x=grille[1]-grille[0]
    
    rho_n=rho_0
    rho_np1=np.zeros(rho_0.shape)

    sol_list=list()
    sol_list.append(rho_0)

    dt=CFL*delta_x/adv
    
    for t in range(nt):
        
        #+++++ Conditions limites +++++
        rho_np1[0]=rho_n[-2]
        rho_np1[-1]=rho_n[1]
        
        #+++++ Interieur du domaine +++++
        
        if stationnaire == 'True':
            
            f_n=rho_n[1:-1]*u_max
            f_np1=rho_n[2:]*u_max
            f_nm1=rho_n[0:-2]*u_max
            
            # Correction d'Harten 
            a_jp1=0.5*(rho_n[2:]+rho_n[1:-1])
            
            for i in range(np.size(a_jp1)):
                epsilon = np.sqrt(np.max([0,(a_jp1[i]-u_max*(1-2*rho_n[i]/rho_max)), (-a_jp1[i]+(u_max*(1-2*rho_n[i]/rho_max)))]))
                if a_jp1[i]<epsilon :
                    a_jp1[i]=0.5*((np.abs(a_jp1[i])**2)/epsilon)+0.5*epsilon
    
            a_jm1=0.5*(rho_n[1:-1]+rho_n[0:-2])
        
        elif stationnaire == 'False':

            f_n=rho_n[1:-1]*u_max*(1-rho_n[1:-1]/rho_max)
            f_np1=rho_n[2:]*u_max*(1-rho_n[2:]/rho_max)
            f_nm1=rho_n[0:-2]*u_max*(1-rho_n[0:-2]/rho_max)
            
            # Correction d'Harten 
            a_jp1=0.5*(rho_n[2:]+rho_n[1:-1])
            
            for i in range(np.size(a_jp1)):
                epsilon = np.sqrt(np.max([0,(a_jp1[i]-u_max*(1-2*rho_n[i]/rho_max)), (-a_jp1[i]+(u_max*(1-2*rho_n[i]/rho_max)))]))
                if a_jp1[i]<epsilon :
                    a_jp1[i]=0.5*((np.abs(a_jp1[i])**2)/epsilon)+0.5*epsilon
           
            a_jm1=0.5*(rho_n[1:-1]+rho_n[0:-2])
            
        else:
            print('Choix du type de calcul invalide')
            break
                        
        fnum_p=(f_np1+f_n)/2-np.abs(a_jp1)/2*(rho_n[2:]-rho_n[1:-1])
        fnum_m=(f_n+f_nm1)/2-np.abs(a_jm1)/2*(rho_n[1:-1]-rho_n[0:-2])
        
        rho_np1[1:-1]=rho_n[1:-1]-dt/delta_x*(fnum_p-fnum_m)
        
        #+++++ Mise a jour de la solution +++++
        rho_n=copy.copy(rho_np1)

        #+++++ Enregistrement de la solution +++++
        sol_list.append(rho_n)

    return np.array(sol_list)


def visualisation(step):
    grille=step['maillage']['grille']
    sol_list=step['calcul']['resultat']
    rho_0=step['calcul']['init']
    
    #----- Figures -----
    figure=plt.figure(figsize=[7.0,5.0])
    plt.xlim([step['maillage']['xmin'],step['maillage']['xmax']])
    plt.ylim([0.2,1.4])
    plt.plot(grille,sol_list[-1],'ro', label ='Solution à t=50')
    plt.plot(grille,rho_0,'k-', label='Solution à t=0')
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')

    plt.savefig(step['visualisation']['Image']['Fichier'])
    plt.show()

#######################
##### SCRIPT PART #####
#######################
#***** Parametres *****
step=dict()

step['maillage']=dict()
step['init']=dict()
step['calcul']=dict()
step['visualisation']=dict()

#----- Maillage -----
step['maillage']['Nx']=101
step['maillage']['xmin']=-2.0
step['maillage']['xmax']=2.0
step['maillage']['grille']=make_grille(step)

#----- Solution initiale -----
step['init']['rho_max']=1.0
step['calcul']['init']=init(step)

#----- Calcul -----
step['calcul']['adv']=1.0
step['calcul']['CFL']=0.5
step['calcul']['nt']=50
step['calcul']['u_max']=1.0
step['calcul']['stationnaire']= 'False' # 'True' ou 'False'

#***** Schema numerique *****
step['calcul']['resultat']=upwind(step)

#***** Sortie graphique *****
step['visualisation']['Image']=dict()
step['visualisation']['Image']['Fichier']='Initiale_finale.png' # jpg, png, pdf, ps, eps et svg.

visualisation(step)
