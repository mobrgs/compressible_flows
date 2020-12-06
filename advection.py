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
def make_grille(paradic):
    xmin=paradic['maillage']['xmin']
    xmax=paradic['maillage']['xmax']
    Nx=paradic['maillage']['Nx']
    return np.linspace(xmin,xmax,Nx)

def init(paradic):
    grille=paradic['maillage']['grille']
    mu=paradic['init']['mu']
    sigma=paradic['init']['sigma']
    return np.exp(-0.5*((grille-mu)/sigma)**2.0)*np.exp(0.5*(mu/sigma)**2.0)

def upwind(paradic):
    w0=paradic['calcul']['init']
    grille=paradic['maillage']['grille']
    nt=paradic['calcul']['nt']
    
    CFL=paradic['calcul']['CFL']
    adv=paradic['calcul']['adv']
    
    delta_x=grille[1]-grille[0]
    
    wn=w0
    wnp1=np.zeros(w0.shape)

    sol_list=list()
    sol_list.append(w0)

    dt=CFL*delta_x/adv
    imax=w0.shape[0]-1 # retourne l'indice de la dernière maille
    
    for t in range(nt):
        
        #+++++ Conditions limites +++++
        wnp1[0]=wn[imax]
        
        #+++++ Interieur du domaine +++++
        wnp1[1:imax]=wn[1:imax]-adv*dt/delta_x*(wn[1:imax]-wn[0:imax-1])
        
        #+++++ Mise a jour de la solution +++++
        wn=copy.copy(wnp1)
        
        print(max(wn)) # Retourne la valeur maximale de la courbe

        #+++++ Enregistrement de la solution +++++
        sol_list.append(wn)

    return np.array(sol_list)

def centred(paradic):
    w0=paradic['calcul']['init']
    grille=paradic['maillage']['grille']
    nt=paradic['calcul']['nt']
    
    CFL=paradic['calcul']['CFL']
    adv=paradic['calcul']['adv']
    theta=paradic['calcul']['theta']
    
    delta_x=grille[1]-grille[0]
    
    wn=w0
    wnp1=np.zeros(w0.shape)

    sol_list=list()
    sol_list.append(w0)

    dt=CFL*delta_x/adv
    imax=w0.shape[0]-1
    
    for t in range(nt):
        #+++++ Conditions limites +++++
        wnp1[0]=wn[imax]
    
        #+++++ Interieur du domaine +++++
        wnp1[1:imax-1]=wn[1:imax-1]+theta/2*(wn[2:imax]-2*wn[1:imax-1]+wn[0:imax-2])-adv*dt/2/delta_x*(wn[2:imax]-wn[0:imax-2])

        #+++++ Conditions limites +++++
        wnp1[imax]=wnp1[0]

        #+++++ Mise a jour de la solution +++++
        wn=copy.copy(wnp1)

        #+++++ Enregistrement de la solution +++++
        sol_list.append(wn)
        
        

    return np.array(sol_list)

def visualisation(paradic):
    grille=paradic['maillage']['grille']
    sol_list=paradic['calcul']['resultat']
    w0=paradic['calcul']['init']
    
    #----- Figures -----
    figure=plt.figure(figsize=[7.0,5.0])
    plt.xlim([paradic['maillage']['xmin'],paradic['maillage']['xmax']])
    plt.ylim([0.0,1.2])
    plt.plot(grille,sol_list[-1],'ro', label ='Solution à t=30')
    plt.plot(grille,w0,'k-', label='Solution à t=0')
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('w')
    plt.title('Theta = 1.0 et CFL = 1.0')

    plt.savefig(paradic['visualisation']['Image']['Fichier'])
    
    #----- Animation -----
    # plot_list=list()
    
    # figure=plt.figure(figsize=[7.0,5.0])
    # plt.xlim([paradic['maillage']['xmin'],paradic['maillage']['xmax']])
    # plt.ylim([0.0,1.2])
    # plt.grid()
    
    # for sol in sol_list:
    #     plot_list.append(plt.plot(grille,sol,'ro',grille,w0,'k-'))

    # animation=anim.ArtistAnimation(figure,plot_list,interval=50)

    #----- Affichage -----
    plt.show()

#######################
##### SCRIPT PART #####
#######################
#***** Parametres *****
paradic=dict()

paradic['maillage']=dict()
paradic['init']=dict()
paradic['calcul']=dict()
paradic['visualisation']=dict()

#----- Maillage -----
paradic['maillage']['Nx']=101
paradic['maillage']['xmin']=-2.0
paradic['maillage']['xmax']=2.0
paradic['maillage']['grille']=make_grille(paradic)

#----- Solution initiale -----
paradic['init']['mu']=0.0
paradic['init']['sigma']=0.4
paradic['calcul']['init']=init(paradic)

#----- Calcul -----
paradic['calcul']['adv']=0.5
paradic['calcul']['CFL']=0.01
paradic['calcul']['nt']=1000
paradic['calcul']['theta']=1.0 # Pour le schema centre uniquement.

#***** Schema numerique *****
paradic['calcul']['resultat']=upwind(paradic)
#paradic['calcul']['resultat']=centred(paradic)

#***** Sortie graphique *****
paradic['visualisation']['Image']=dict()
paradic['visualisation']['Image']['Fichier']='Initiale_finale.png' # jpg, png, pdf, ps, eps et svg.

visualisation(paradic)
