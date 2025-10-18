#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 21:46:51 2022

@author: rdempsey
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:30:49 2022

@author: rdempsey
"""

import numpy as np
import os
import itertools, datetime
from math import isclose

"some handy functions for translating between units"

h_bar=6.5821220*(10**(-16))
q_e=1.60217733*(10**(-19))

def beta_to_temp(beta):
    return(1/(beta*8.617385*(10**(-5))))

def temp_to_beta(temp):
    return(1/(temp*8.617385*(10**(-5))))

def picosecs(time):
    return(time*6.5821220*10**(-4))

def timeticks(picoseconds):
    return(picoseconds/(6.5821220*10**(-4)))

def microampere(current):
    return((10**(6))*current*q_e/h_bar)

def microampere_squared(currentsquared):
    return(currentsquared*((q_e/h_bar)**2)*(10**(12)))

"translates site and band number into the apppriate basis number"

def basis(L,n):
    d=2*(2*L+1)
    v=np.zeros(d)
    v[n]=1
    return v

def index(site,band,L):
    return(site+L+band*(2*L+1))

"defines some auxillary observables useful for defining the Hamiltonian"
    
def count_operator(site,band,L):
    return np.outer(basis(L,index(site,band,L)), basis(L,index(site,band,L)))

def hop_left(site,band,L):
    return np.outer(basis(L,index(site-1,band,L)), basis(L,index(site,band,L)))

def hop_right(site,band,L):
    return np.outer(basis(L,index(site+1,band,L)), basis(L,index(site,band,L)))

def hop_band(site,L):
    return np.outer(basis(L,index(site,0,L)), basis(L,index(site,1,L)))+ np.outer(basis(L,index(site,1,L)), basis(L,index(site,0,L)))

"gives back the discrete Laplacian"

def D_Laplacian(s,ran,L):
    return (s*(2*sum([count_operator(site,1,L) for site in range(ran[0],ran[1])])
            -sum([hop_left(site,1,L) for site in range(ran[0]+1,ran[1])])
            -sum([hop_right(site,1,L) for site in range(ran[0],ran[1]-1)])))

"Hamilitonian of the system without the voltage applied"

def H_0(epsilon, kappa, theta, mu, gamma, l,L):
    
    "define the kinetic energy part"
    
    H_K_pro=lambda epsilon,l,L : D_Laplacian(epsilon,(-l,l+1),L)

    H_K_res_l=lambda kappa,l,L : D_Laplacian(kappa,(-L,-l),L)
    H_K_res_r=lambda kappa,l,L : D_Laplacian(kappa,(l+1,L+1),L)

    H_K_res_pro=lambda theta,l,L: -theta*(hop_left(-l,1,L)+hop_right(-l-1,1,L)+hop_left(l+1,1,L)+hop_right(l,1,L))

    "define the chemical potential"
    H_chem=lambda mu,l,L: (mu[0]*sum([count_operator(site,0,L) for site in range(-L,L+1)])
                                   +mu[1]*sum([count_operator(site,1,L) for site in range(-l,l+1)])
                                   +mu[2]*sum([count_operator(site,1,L) for site in range(-L,-l)])
                                   +mu[3]*sum([count_operator(site,1,L) for site in range(l+1,L+1)]))

    H_band_hop=lambda gamma,l,L: gamma*sum([hop_band(site,L) for site in range(-l,l+1)])

    return (H_K_pro(epsilon,l,L)+H_K_res_l(kappa,l,L)+H_K_res_r(kappa,l,L)+H_K_res_pro(theta,l,L)
                                                    -H_chem(mu,l,L)-H_band_hop(gamma,l,L))  
  

"Hamilitonian of the system without the voltage applied"

def H_eta(eta,l,L) :
    H_eta=-eta*sum([count_operator(site,1,L) for site in range(-L,-l)])
    H_eta+=eta*sum([(site/l)*count_operator(site,1,L)+(site/l)*count_operator(site,0,L) for site in range(-l,l+1)])
    H_eta=H_eta+eta*sum([count_operator(site,1,L) for site in range(l+1,L+1)])
    return H_eta
    
"Helper functions required to caluculate observables/apply the Hamiltonian"

def apply_rho(E,V,beta,psi):
    W=np.matmul(np.conj(V),psi)*(1/(1+np.exp(E*beta)))
    D=np.matmul(W, V)
    return D;

def apply_Uplus(t,E,V,psi):
    W=np.matmul(np.conj(V),psi)*(np.exp(1j*t*E))
    D=np.matmul(W, V)
    return D;

"function calculating the current evolution"

def current_evolve(t,eigeninfo_0,eigeninfo_eta,l,L,beta):
    E_nu=eigeninfo_eta[0];
    V_nu=eigeninfo_eta[1].T;
    E_0=eigeninfo_0[0];
    V_0=eigeninfo_0[1].T;
    C=0
    for i in range(-l,l+1):
        ex=np.zeros(2*(2*L+1));
        ex1=np.zeros(2*(2*L+1));
        ex[i+3*int(L)+1]=1;
        ex1[i+3*int(L)+2]=1;
        C_t=2*np.imag((np.vdot(apply_Uplus(t,E_nu, V_nu,ex),apply_rho(E_0,V_0,beta,
                                                  apply_Uplus(t,E_nu,V_nu,ex1)))))
        C=C+C_t; 
        R=np.real(C/(2*l))
        
    return R

"Finds the time evolution of the squared current"

def squared_current_alt(t,eigeninfo_0,eigeninfo_eta,current_squared,l,L,beta):
    E_nu=eigeninfo_eta[0];
    V_nu=eigeninfo_eta[1].T;
    E_0=eigeninfo_0[0];
    V_0=eigeninfo_0[1].T;

    C2=0
    for item in current_squared:
        x=basis(L,int(item[1]))
        y=basis(L,int(item[2]))
        u=basis(L,int(item[3]))
        v=basis(L,int(item[4]))
        term1=np.vdot(apply_Uplus(t,E_nu,V_nu,x),apply_rho(E_0,V_0,beta,
                                              apply_Uplus(t,E_nu,V_nu,y)))
        term2=np.vdot(apply_Uplus(t,E_nu,V_nu,u),apply_rho(E_0,V_0,beta,
                                              apply_Uplus(t,E_nu,V_nu,v)))
        term3=np.vdot(apply_Uplus(t,E_nu,V_nu,u),apply_rho(E_0,V_0,beta,
                                              apply_Uplus(t,E_nu,V_nu,y)))
        term4=np.vdot(apply_Uplus(t,E_nu,V_nu,x),apply_rho(E_0,V_0,beta,
                                              apply_Uplus(t,E_nu,V_nu,v)))
        evolved_y=apply_Uplus(t,E_nu,V_nu,y)
        evolved_u=apply_Uplus(t,E_nu,V_nu,u)
        increment=item[0]*(term1*term2-term3*term4+np.vdot(evolved_y,evolved_u)*term4)
        C2=C2+increment

    return C2/((2*l)**2)

"Auxilary functions used to find the form of the sqaured current density"

def sym_current_observable(site,L):
    n=index(site,1,L)
    return([-1j,n,n-1],[1j,n-1,n])

def sym_current_density_squared(l,L):
    s_c_d=[x for site in range(-l+1,l) for x in sym_current_observable(site,L)]
    return([[(x[0]*y[0])]+x[1:]+y[1:] for x in s_c_d for y in s_c_d])



"function calculating the charge transport at a site"

def charge_transport(t,eigeninfo_0,eigeninfo_eta,site,L,beta):
    E_nu=eigeninfo_eta[0];
    V_nu=eigeninfo_eta[1].T;
    E_0=eigeninfo_0[0];
    V_0=eigeninfo_0[1].T;
    C_t=[]
    for ind in [index(site,0),index(site,1)]:

        ex=np.zeros(2*(2*L+1));
        ex1=np.zeros(2*(2*L+1));
        ex[ind]=1;
        ex1[ind+1]=1;
        C_t.append(2*np.imag((np.vdot(apply_Uplus(t,E_nu, V_nu,ex),apply_rho(E_0,V_0,beta,
                                                  apply_Uplus(t,E_nu,V_nu,ex1))))))

    return C_t

"function calculating the particle density"

def particle_density(t,eigeninfo_0,eigeninfo_eta,site,L,beta):
    E_nu=eigeninfo_eta[0];
    V_nu=eigeninfo_eta[1].T;
    E_0=eigeninfo_0[0];
    V_0=eigeninfo_0[1].T;
    C_t=[]
    for ind in [index(site,0),index(site,1)]:
        ex=np.zeros(2*(2*L+1));
        ex1=np.zeros(2*(2*L+1));
        ex[ind]=1;
        ex1[ind]=1;
        C_t.append(np.vdot(apply_Uplus(t,E_nu, V_nu,ex),apply_rho(E_0,V_0,beta,
                                                  apply_Uplus(t,E_nu,V_nu,ex1))))
    return C_t


"fuction to write a textfile containing the paramteers of the calculation"

def gen_metadata(filename, dictionary):
    
    with open(filename+'.txt', 'w') as f:
        for (k,v) in dictionary.items():
            f.write('%s=%s\n' % (k, v))
        f.close()

"if the user wishes the program to find the stationary regime, then AutoSR is called"
        
def autosr(eigeninfo_0,eigeninfo_eta, epsilon, kappa, theta, mu, gamma,eta,l,L,beta,time=0):
    limit=(10**(-3))
    "choose reasonable times"
    times=np.linspace(time,time+100,20)
    I=np.array([microampere(current_evolve(t,eigeninfo_0,eigeninfo_eta,l,L,beta)) for t in times])   
       
    n_L=L+50
    h_0=H_0(epsilon, kappa, theta, mu, gamma, l,n_L)
    n_eigeninfo_0=np.linalg.eig(h_0)
    h_eta=h_0+H_eta(eta,l,n_L)
    n_eigeninfo_eta=np.linalg.eig(h_eta)
    testI=np.array([microampere(current_evolve(t,n_eigeninfo_0,n_eigeninfo_eta,l,n_L,beta)) for t in times])
    I2=testI[1:]
    I1=testI[:-1] 
    
    "while increasing L doesn't do anything increase the time and check"
    while  not np.mean(np.abs(np.abs(I2)-np.abs(I1)))<limit:
        
        if not all([isclose(I[i],testI[i],rel_tol=0.0001) for i in range(10)]):
                print(L)
                (I,L,eigeninfo_0,eigeninfo_eta)=(testI,n_L,n_eigeninfo_0,n_eigeninfo_eta)
                n_L=L+50
                h_0=H_0(epsilon, kappa, theta, mu, gamma, l,n_L)
                n_eigeninfo_0=np.linalg.eig(h_0)
                h_eta=h_0+H_eta(eta,l,n_L)
                n_eigeninfo_eta=np.linalg.eig(h_eta)
                testI=np.array([microampere(current_evolve(t,n_eigeninfo_0,n_eigeninfo_eta,l,n_L,beta)) for t in times])
        print(times)
        times=times+50
        I=np.array([microampere(current_evolve(t,eigeninfo_0,eigeninfo_eta,l,L,beta)) for t in times])
        testI=np.array([microampere(current_evolve(t,n_eigeninfo_0,n_eigeninfo_eta,l,n_L,beta)) for t in times])
        I2=testI[1:]
        I1=testI[:-1]
        
    return testI, n_L, times, n_eigeninfo_0, n_eigeninfo_eta 
                
        
"""Takes in the parameters of the system and what to calculate (Current, Current Squared, Quantum Fluctations, Particle Density). AutoSR
asks if we wish the program to automatically register the stationary regime"""

def Evolution(Ls=100, ls=2, epsilons=0.65, kappas=0.85, thetas=0.75, mu_0s=0, mu_12s=0.5, mu_11s=0.5, gammas=0, etas=0.1, temps=300, C=False, CT=False, CS=False, QF=False, PD=False, AutoSR=False, dirname=None):

    
    "create a directory where we will store the results"
    
    if not dirname:
        today=str(datetime.date.today()).replace('-','')
        dirname=today+'_quasifree_calcs'
    if os.path.exists(dirname):
        dirname=input('directory already exists; please input a name for the new directory\n')
    
    os.mkdir(dirname)
    
    "take in list of sites at which to calculate PD or CT"
    if PD==True:
        pd_sites=input('please input a list of the sites at which the particle density is to be calculated')

    if CT==True:
        ct_sites=input('please input a list of the sites at which the particle density is to be calculated')      
    "convert units; convert units"
    
    "if the user wants to input the times manually they should supply a numpy file with times"
    
    if AutoSR==False:
        fn=input('please input the filename of a .npy file with desired times\n')
        times=timeticks(np.load(fn))
    else:
        
        "otherwise initialize lists and dictionaries to hold average of observables, the related times and Ls"

        (aLs,atimes,meanIs,meanI2s)=([],[],[],[])
        if PD==True:
            dict_pd={site: [] for site in pd_sites}
        if CT==True:
            dict_ct={site: [] for site in ct_sites} 
        
    list_parameters=[]
    
    for x in [Ls, ls, epsilons, kappas, thetas, mu_0s, mu_12s, mu_11s, mu_11s, gammas]:
        if type(x)==list:
            list_parameters.append(list(x))
        else:
            list_parameters.append([x])
            
    if not type(etas)==list:
        etas=[etas]
    
    if not type(temps)==list:
        temps=[temps]
        
    betas=[]
    for x in temps:
        betas.append(temp_to_beta(x))

    "form a list of the combinations of all 'fundamental' quantities"

    c0=list(itertools.product(*list_parameters))

    "starting counter to count files/generate unique filenames"
    
    i=0
    
    "The desired calculations are done for all combinations of parameters."
    print(len(etas))
    for params in list(c0):

        (L,l,epsilon, kappa, theta, mu_0, mu_12, mu_11, mu_11, gamma)=params
        mu=[mu_0, mu_12, mu_11, mu_11]
        h_0=H_0(epsilon, kappa, theta, mu, gamma, l,L)
        eigeninfo_0=np.linalg.eig(h_0)
        print(type(h_0))

        for eta in etas:
            h_eta=h_0+H_eta(eta,l,L)
            eigeninfo_eta=np.linalg.eig(h_eta)

            for beta in betas:
                
                i=i+1

                print('i is ',i)
                "generate a metadata file giving the params involved in the calculations"
                filename='parameters_'+str(i)
                filename='./'+dirname+'/'+filename
                param_titles=['L', 'l', 'epsilon', 'kappa', 'theta','mu1','mu2','mu3','mu4', 'gamma','eta','beta','AutoSR']
                fparams=list(params)+[eta,beta,AutoSR]
                dict_metadata={param_titles[i]:fparams[i] for i in range(len(fparams))}
                gen_metadata(filename, dict_metadata)
                if AutoSR==True:
                        (Is,L,times,eigeninfo_0,eigeninfo_eta)=autosr(eigeninfo_0,eigeninfo_eta, epsilon, kappa, theta, mu, gamma,eta,l,L,beta)
                        meanIs.append(np.mean(Is))
                        np.save(dirname+'/current_%s.npy' % i,np.array([picosecs(times),Is]))
                        aLs.append(L)
                        atimes.append(picosecs(times))
                        h_0=H_0(epsilon, kappa, theta, mu, gamma, l,L)
                        h_eta=h_0+H_eta(eta,l,L)
                        with open(dirname+'/SR_parameters.txt', 'w') as f:
                            f.write('These are the parameters for which a stationary regime was found\n')
                            f.write('L=%s\n' % aLs)
                            f.write('times=%s\n' % atimes)
                            f.close()
                            
                        
                if C==True and AutoSR==False: 
                        I=[microampere(current_evolve(t,eigeninfo_0,eigeninfo_eta,l,L,beta)) for t in times]
                        np.save(dirname+'/current_%s.npy' % i,np.array([picosecs(times),I]))

                                    
                if CS==True:
                    current_squared=sym_current_density_squared(l,L)
                    if AutoSR==True:
                        I2=[microampere_squared(squared_current_alt(t,eigeninfo_0,eigeninfo_eta,current_squared,l,L,beta)) for t in times]
                        meanI2s.append(np.mean(I2))
                        np.save(dirname+'/current_squared_%s.npy' % i,np.array([picosecs(times),I2]))
                    else: 
                        I2=[microampere_squared(squared_current_alt(t,eigeninfo_0,eigeninfo_eta,current_squared,l,L,beta)) for t in times]
                        np.save(dirname+'/current_squared_%s.npy' % i,np.array([picosecs(times),I2]))


                if QF==True:
                    if AutoSR==True:
                        if not CS:
                            current_squared=sym_current_density_squared(l,L)
                            I2=[microampere_squared(squared_current_alt(t,eigeninfo_0,eigeninfo_eta,current_squared,l,L,beta)) for t in times]
                            np.save(dirname+'/current_squared_%s.npy' % i,np.array([picosecs(times),I2]))
                            meanI2s.append(np.mean(I2))

                    else: 
                        if not CS: 
                            current_squared=sym_current_density_squared(l,L)
                            I2=[microampere_squared(squared_current_alt(t,eigeninfo_0,eigeninfo_eta,current_squared,l,L,beta)) for t in times]
                        if not C:
                            I=[microampere(current_evolve(t,eigeninfo_0,eigeninfo_eta,l,L,beta)) for t in times]
                            
                        variance=np.array(I2)-np.array(I)**2
                        filename=dirname+'/variance_%s.npy' % i
                        np.save(filename,np.array([picosecs(times),variance]))
                        
                if PD==True:
                    sites=input('please input a list of the sites at which the particle density is to be calculated')
                    if AutoSR==True:
                        for site in pd_sites:
                            dict_pd[site].append(np.mean(np.array([particle_density(t,eigeninfo_0,eigeninfo_eta,site,L,beta) for t in times])))
                    else: 
                        for site in pd_sites:
                            PD_site=np.array([particle_density(t,eigeninfo_0,eigeninfo_eta,site,L,beta) for t in times])
                            filename=dirname+'/PD_site_%s_%s.npy' % (site, i)
                            np.save(filename,np.array(picosecs(times),PD_site))
                            
                if CT==True:
                    sites=input('please input a list of the sites at which the charge transport is to be calculated')
                    if AutoSR==True:
                        for site in ct_sites:
                            dict_ct[site].append(np.mean(np.array([charge_transport(t,eigeninfo_0,eigeninfo_eta,site,L,beta) for t in times])))
                    else: 
                        for site in ct_sites:
                            CT_site=[charge_transport(t,eigeninfo_0,eigeninfo_eta,site,L,beta) for t in times]
                            filename=dirname+'/CT_site_%s_%s.npy' % (site, i)
                            np.save(filename,np.array(picosecs(times),CT_site))
                            
                    
          
                "If AutoSR was chosen write all the information to the chosen directory"
                
                if AutoSR:
            
                    if C: np.save(dirname+'/current.npy',meanIs)
                    if CS: np.save(dirname+'/current_squared.npy',meanI2s)
                    if QF: np.save(dirname+'/variance.npy',meanI2s-np.array(meanIs)*np.array(meanIs))
                    if CT: 
                        for site in ct_sites:
                            np.save(dirname+'/charge_transport_site_%s.npy' % sites,dict_ct[site])
                    if PD: 
                        for site in pd_sites:
                            np.save(dirname+'/particle_density_site_%s.npy' % sites,dict_ct[site])

