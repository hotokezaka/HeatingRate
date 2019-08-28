from numpy import pi
import numpy as np
import pandas as pd
import bateman as bt
import thermalization as th
from scipy import optimize
day = 86400.
eV = 1.60218e-12
MeV = 1.0e6*eV
mu = 1.66054e-24
c = 2.99792458e10
def calc_heating_rate(Mej,vej, Amin,Amax,ffraction,kappa_effs,alpha_max,alpha_min,n):

    t_initial = 0.01*day
    t_final = 1000.*day
    delta_t = 0.3#0.3

    Nth = 40#40 for default

    
    Amax_beta = 209
    
    Xtot = 0.0

    fraction = np.zeros(300)
    for i in range(0,len(ffraction)):
        A = ffraction[1][i]
        fraction[A] = float(A)*ffraction[2][i]
    tmpA = 0.0
    for A in range(Amin,Amax+1):
        Xtot+=fraction[A]
        tmpA += float(A)*fraction[A]

    Aave = tmpA/Xtot



###    




    total_heats = []
    total_gammas = []
    total_elects = []
    total_elect_ths = []
    total_gamma_ths = []
    heating_functions = []
    ts = []

    t = t_initial
    while t<t_final:
        total_heats.append(0.)
        total_gammas.append(0.)
        total_elects.append(0.)
        total_elect_ths.append(0.)
        total_gamma_ths.append(0.)
        heating_functions.append(0.)
        ts.append(t)
        t*= 1. + delta_t

    print 'total time step = ', len(total_heats)
    for A in range(Amin,min(Amax,Amax_beta)+1): 
        each_heats = np.zeros(len(ts))
        each_gammas = np.zeros(len(ts))
        each_gamma_ths = np.zeros(len(ts))
        each_elects = np.zeros(len(ts))
        each_elect_ths = np.zeros(len(ts))
        Xfraction = fraction[A]/Xtot
    
        filename = 'input_files/table_beta/'+str(A)+'.txt'
        filename2 = 'heat'+str(A)+'.dat'
    
#A, Z, Q[MeV], Egamma[MeV], Eelec[MeV], Eneutrino[MeV], tau[s]
        fchain = pd.read_csv(filename,delim_whitespace=True,header=None)
#    print 'length of the chain',A,len(fchain)

        tmp = []
        N = len(fchain)
        for i in range(0,N):
            tmp.append(1.0/fchain[6][i])
        lambdas = np.array(tmp)
####determine the thermalization time in units of day for each element
        tes = np.zeros(N)
        total_numb = np.zeros(N)
        for i in range(0,N):
            Z = fchain[1][i]
            Egamma = fchain[3][i]
            Eele = fchain[4][i]
       
            if(Eele>0.):
                tes[i] = th.calc_thermalization_time(Eele,Mej,vej,Aave,alpha_max,alpha_min,n)


    #ts = []
        tmp_numb = np.zeros(N)
        number10s = []
        number11s = []
        number12s = []
        number13s = []
        number14s = []
        number15s = []

        heats = []
        gammas = []
        elects = []

        lambda_sort = np.zeros((N,N))
        coeffs = np.ones(N)
        xs = np.zeros((N,N))

        lambda_sort = np.zeros((N,N))

        for i in range(0,N):
            tmp0 = lambdas[:i+1]
            tmp = np.sort(tmp0)
    
            lambda_sort[i][:i+1] = tmp[::-1]

        for j in range(1,N):
            for i in range(0,j):
                coeffs[j] *= lambda_sort[j-1][i]
   # print j,coeffs[j]

        for j in range(1,N):
            for i in range(0,j):
                xs[j][i] = (lambda_sort[j][j-1-i]-lambda_sort[j][j])
    
    
        for k in range(0,len(ts)):
   # while t<t_final:
            t = ts[k]
        
            for i in range(0,N):
                coeff = coeffs[i]*np.power(t,i)
        
                if(i==6):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_6(xs[6][0]*t,xs[6][1]*t,xs[6][2]*t,xs[6][3]*t,xs[6][4]*t,xs[6][5]*t)

                elif(i==5):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_5(xs[5][0]*t,xs[5][1]*t,xs[5][2]*t,xs[5][3]*t,xs[5][4]*t)
          
                elif(i==4):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_4(xs[4][0]*t,xs[4][1]*t,xs[4][2]*t,xs[4][3]*t)
          #  print i,coeff
                elif(i==3):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_3(xs[3][0]*t,xs[3][1]*t,xs[3][2]*t)
          #  print i,xs[3][0],xs[3][1],xs[3][2]
                elif(i==2):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_2(xs[2][0]*t,xs[2][1]*t)
          #  print i,xs[2][0],xs[2][1]
                elif(i==1):
                    tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_1(xs[1][0]*t)
          #  print i,xs[1][0]
                elif(i==0):
                    tmp_numb[i] = np.exp(-t*lambda_sort[i][i])
                else:
                    print 'chain is too long'
           # print i,lambda_sort[i][i]
            
        
#        number10s.append(tmp_numb[0])
#        number11s.append(tmp_numb[1])
#        number12s.append(tmp_numb[2])
#        number13s.append(tmp_numb[3])
    
            heat = 0.0
            gam = 0.0
            ele = 0.0
            ele_th = 0.0
            gam_th = 0.0
    
            for i in range(0,N):
                
                Eele = fchain[4][i]
                if(t > 0.003*tes[i]):
                    if(Eele > 0.):
                        tau1 = t/tes[i]
                        if(tau1<2.):
                            tau0 = 0.03*tau1#0.05*tau1
                            root = optimize.newton(th.calc_zero_energy, tau0, args=(tau1,Eele,))
                            tau0 = root
                        else:
                            tau0 = 0.4#0.05*tau1
                #else:
                #    tau0 = 0.01*tau1
                #tau0 = 1.
                
                #print "tau0: ", tau1,tau0
                        delta_t = (tau1-tau0)*tes[i]
            #fth = fchain[6][i]*tes[i]*tes[i]*(np.exp(delta_t/fchain[6][i])-1.0)*np.power(t,-3.)
           # print t,delta_t,delta_t/t,fchain[6][i]*tes[i]*tes[i]*(np.exp(delta_t/fchain[6][i])-1.0)*np.power(t,-3.)
                        t_th = tau0*tes[i]
                        tmp_n_t_th = 1./float(Nth-1)
                        dt_th = np.power(tau1/tau0,tmp_n_t_th)-1.
                        total_numb[i] = 0.0
                        for j in range(0,Nth):
                            coeff = coeffs[i]*np.power(t_th,i)
                            tau = t_th/tes[i]
                    #e_delay= calc_e_tau_tau0(tau,tau1)
                            e_delay = th.epsilon_tau(tau,tau1,Eele)

                            if(i==6):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_6(xs[6][0]*t_th,xs[6][1]*t_th,xs[6][2]*t_th,xs[6][3]*t_th,xs[6][4]*t_th,xs[6][5]*t_th)

                            elif(i==5):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_5(xs[5][0]*t_th,xs[5][1]*t_th,xs[5][2]*t_th,xs[5][3]*t_th,xs[5][4]*t_th)
          
                            elif(i==4):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_4(xs[4][0]*t_th,xs[4][1]*t_th,xs[4][2]*t_th,xs[4][3]*t_th)
                            elif(i==3):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_3(xs[3][0]*t_th,xs[3][1]*t_th,xs[3][2]*t_th)
                            elif(i==2):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_2(xs[2][0]*t_th,xs[2][1]*t_th)
                            elif(i==1):
                                tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_1(xs[1][0]*t_th)
                            elif(i==0):
                                tmp_numb[i] = e_delay*np.exp(-t_th*lambda_sort[i][i])
                            

                            total_numb[i] += tmp_numb[i]*dt_th*t_th
                                                                                  
                            t_th = (1.+dt_th)*t_th
                    if(Egamma>0.):
                        Z = fchain[1][i]
                        kappa_eff = kappa_effs[A][Z]
                        fth_gamma = th.calc_gamma_deposition(kappa_eff,t,Mej,vej,alpha_min,alpha_max,n)
                    else:
                        fth_gamma = 0.
                
                    heat += Xfraction*MeV*tmp_numb[i]*fchain[2][i]*lambdas[i]/(mu*float(A))
                    gam += Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A))
                    gam_th += fth_gamma*Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A))
                    ele += Xfraction*MeV*tmp_numb[i]*fchain[4][i]*lambdas[i]/(mu*float(A))
                    ele_th += np.power(tes[i],2.)*np.power(t,-3.)*Xfraction*MeV*total_numb[i]*fchain[4][i]*lambdas[i]/(mu*float(A))           
                else:
                        

           
                    if(Egamma>0.):
                        Z = fchain[1][i]
                        kappa_eff = kappa_effs[A][Z]
                        fth_gamma = th.calc_gamma_deposition(kappa_eff,t,Mej,vej,alpha_min,alpha_max,n)
                    else:
                        fth_gamma = 0.
                
                    heat += Xfraction*MeV*tmp_numb[i]*fchain[2][i]*lambdas[i]/(mu*float(A))
                    gam += Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A))
                    gam_th += fth_gamma*Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A))
                    ele += Xfraction*MeV*tmp_numb[i]*fchain[4][i]*lambdas[i]/(mu*float(A))
                    ele_th += Xfraction*MeV*tmp_numb[i]*fchain[4][i]*lambdas[i]/(mu*float(A))



            total_heats[k] += heat
            total_gammas[k] += gam

            total_elects[k] += ele
            total_elect_ths[k] += ele_th
            total_gamma_ths[k] += gam_th
        


            each_heats[k] += heat
            each_gammas[k] += gam
            each_elects[k] += ele
            each_elect_ths[k] += ele_th
            each_gamma_ths[k] += gam_th
#        print A, Xfraction
    data = {'t':np.multiply(ts,1./day),'total':total_heats,'gamma':total_gammas, 'electron':total_elects, 'gamma_th':total_gamma_ths,'electron_th':total_elect_ths}
    return data        
       # t *= 1.0 + delta_t
    print 'end'


def calc_heating_rate_alpha(Mej,vej, Amin,Amax,ffraction,Yas,kappa_effs,alpha_max,alpha_min,n):
    t_initial = 0.01*day
    t_final = 3000.*day
    delta_t = 0.3

    Nth = 70

    Xtot = 0.

    fraction = np.zeros(300)
    for i in range(0,len(ffraction)):
        A = ffraction[1][i]
        fraction[A] = float(A)*ffraction[2][i]
    tmpA = 0.0
    for A in range(Amin,Amax+1):
        Xtot+=fraction[A]
        tmpA += float(A)*fraction[A]

    Aave = tmpA/Xtot


    total_heats = []
    total_gammas = []
    total_elects = []
    total_alphas = []
    total_elect_ths = []
    total_gamma_ths = []
    total_alpha_ths = []
    ts = []

    t = t_initial
    while t<t_final:
        total_heats.append(0.)
        total_gammas.append(0.)
        total_gamma_ths.append(0.)
        total_elects.append(0.)
        total_elect_ths.append(0.)
        total_alphas.append(0.)
        total_alpha_ths.append(0.)
        ts.append(t)
        t*= 1. + delta_t


##start alpha chain
#A, Z, Q[MeV], Egamma[MeV], Eelec[MeV], Eneutrino[MeV], tau[s]
#    Yas = np.zeros(239)
#    for Aa_tmp in range(210,238):
#        if(Aa_tmp == 222):
#            Yas[Aa_tmp] = 4.0e-5
#        if(Aa_tmp == 223):
#            Yas[Aa_tmp] = 2.7e-5
#        if(Aa_tmp == 224):
#            Yas[Aa_tmp] =  4.1e-5
#        if(Aa_tmp == 225):
#            Yas[Aa_tmp] = 2.7e-5
        
        
    for Aa in range(210,238):
        if(Yas[Aa]!=0.0):
            filename = 'input_files/table_alpha/short_chain_A'+str(Aa)+'.dat'
            fchain = pd.read_csv(filename,delim_whitespace=True,header=None)
    ################
            number_of_br = 0
            for i in range(0,len(fchain)-1):
                if(fchain[2][i]>0.99 or fchain[3][i]>0.99):
                    number_of_br += 0
                else:
                    number_of_br += 1
            number_of_chain = np.power(2,number_of_br)
            print Aa," number of chains:", number_of_chain

#identify all chains starting with A, Z
#the number of chain is the number of branching decays +1
#A=210 is special I take its beta branching to be 1

            number_of_br = 0
            for i in range(0,len(fchain)-1):
                if(fchain[2][i]>0.99 or fchain[3][i]>0.99):
                    number_of_br += 0
                else:
                    number_of_br += 1
            number_of_chain = np.power(2,number_of_br)
        #print "branching points:", number_of_chain





######make chains
            N = len(fchain)
            lambdas = np.zeros((number_of_chain,N))
            Q_betas = np.zeros((number_of_chain,N))
            Q_alphas = np.zeros((number_of_chain,N))
            Q_gammas = np.zeros((number_of_chain,N))
            Q_total_alphas = np.zeros((number_of_chain,N))
            Q_total_betas = np.zeros((number_of_chain,N))
            br_betas = np.zeros((number_of_chain,N))
            br_alphas = np.zeros((number_of_chain,N))
            branching_ratio = np.ones((number_of_chain,N))
            fraction_of_chain = np.zeros((number_of_chain,N))
            length_of_chains = []
            member_of_chain = np.zeros((number_of_chain,N))

            for k in range(0,number_of_chain):
                i = 0
                A = fchain[0][0]
                Z = fchain[1][0]
                count = 0
                count2 = 0
                count3 = 0
                length = 0
                while(fchain[4][i]>0.0):
                    member_of_chain[k][length] = i
                    br_betas[k][length] = fchain[3][i]
                    br_alphas[k][length] = fchain[2][i]
                    Q_total_alphas[k][length] = fchain[5][i]
                    Q_total_betas[k][length] = fchain[6][i]
                    if(fchain[3][i]!=0.):
                        Q_betas[k][length] = fchain[9][i]/fchain[3][i]
                        Q_gammas[k][length] = fchain[8][i]/fchain[3][i]
                    else:
                        Q_betas[k][length] = 0.
                        Q_gammas[k][length] = 0.
                    if(fchain[2][i]!=0.):
                        Q_alphas[k][length] = fchain[7][i]/fchain[2][i]
                    else:
                        Q_alphas[k][length] = 0.

        
                    lambdas[k][length] =np.log(2.)/fchain[4][i]
                    if(number_of_chain==1):
                        fraction_of_chain[k][length] = 1.0
                        if(fchain[2][i]>0.99):
                            A = A -4
                            Z = Z-2
                        #print "alpha", k, A+4, Z+2, "->", A, Z
                        elif(fchain[3][i]>0.99):
                            Z = Z +1
                        #print "beta", k, A, Z-1, "->", A, Z
           
        
                    if(number_of_chain==2):
                        if(fchain[2][i]>0.99):
                            A = A -4
                            Z = Z-2
                            if(count==0):
                                fraction_of_chain[k][length] = 0.5
                    
                            else:
                                fraction_of_chain[k][length] = 1.0
                       # print "alpha", k, A+4, Z+2, "->", A, Z
                        elif(fchain[3][i]>0.99):
                            Z = Z +1
                        #print "beta", k, A, Z-1, "->", A, Z
                            if(count==0):
                                fraction_of_chain[k][length] = 0.5
                    
                            else:
                                fraction_of_chain[k][length] = 1.0
                        elif(k==0):
                        #print "First Branch (alpha)", k,  A+4, Z+2, "->", A, Z
                            A = A -4
                            Z = Z-2
                            branching_ratio[k][i] = fchain[2][i]
                            count +=1
                            fraction_of_chain[k][length] = 1.0
                        elif(k==1):
                        #print "First Branch (beta)", k, A, Z-1, "->", A, Z
                            Z = Z+1
                            count +=1
                            branching_ratio[k][i] = fchain[3][i]
                            fraction_of_chain[k][length] = 1.0
        
                    if(number_of_chain==4):
                        if(fchain[2][i]>0.99):
                            A = A -4
                            Z = Z-2
                       # print "alpha", k, A+4, Z+2, "->", A, Z
                            if(count2==0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                        elif(fchain[3][i]>0.99):
                            Z = Z +1
                        #print "beta", k, A, Z-1, "->", A, Z
                            if(count2==0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                        elif(k==0):
                #alpha-alpha
                        #print "First Branch (alpha)", k,  A+4, Z+2, "->", A, Z
                            A = A -4
                            Z = Z-2
                
                            if(count2 == 0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                            count2 +=1##
                            branching_ratio[k][i] =fchain[2][i]
                        elif(k==1):
                #alpha-beta
                            if(count2 == 0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                            count2 +=1
                
                            if(count == 1):
                            #print "First Branch (alpha)", k,  A+4, Z+2, "->", A, Z
                                A = A -4
                                Z = Z-2
                                branching_ratio[k][i] = fchain[2][i]
                                count += 1
                            else:
                            #print "Second Branch (beta)", k, A, Z-1, "->", A, Z
                                Z = Z+1
                                branching_ratio[k][i] = fchain[3][i]
                        elif(k==2):
                            if(count2 == 0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                            count2 +=1
                #beta-alpha
                            if(count == 0):
                            #print "First Branch (beta)", k,  A+4, Z+2, "->", A, Z
                                Z = Z+1                    
                                count += 1
                                branching_ratio[k][i] = fchain[3][i]
                    
                            else:
                            #print "Second Branch (alpha)", k, A, Z-1, "->", A, Z
                                A = A -4
                                Z = Z-2
                                branching_ratio[k][i] = fchain[2][i]
                    
                        else:
                #beta-beta
                        #print "Branch (beta)", k, A, Z-1, "->", A, Z
                            Z = Z+1
                            branching_ratio[k][i] = fchain[3][i]
                            if(count2 == 0):
                                fraction_of_chain[k][length] = 0.25
                            elif(count2 == 1):
                                fraction_of_chain[k][length] = 0.5
                            else:
                                fraction_of_chain[k][length] = 1.0
                            count2 +=1
                
         
        
                    j = 0
        
                    while(A!=fchain[0][j] or Z!=fchain[1][j]):
                        j += 1
                    i = j
                    length +=1
                length_of_chains.append(length)  
    
    
    ####end make chain############
    
            NN = length_of_chains[0] 
            branch_chains = np.zeros((number_of_chain,NN))
            for k in range(0,number_of_chain):
                branch = 1.0
                for i in range(0,NN):
                    branch *= branching_ratio[k][i]
                    branch_chains[k][i] = branch
            print 'total time step = ', len(total_heats)
            for kk in range(0,number_of_chain): 
                each_heats = np.zeros(len(ts))
                each_gammas = np.zeros(len(ts))
                each_gamma_ths = np.zeros(len(ts))
                each_elects = np.zeros(len(ts))
                each_elect_ths = np.zeros(len(ts))
                each_alphas = np.zeros(len(ts))
                each_alpha_ths = np.zeros(len(ts))
                Xfraction = float(Aa)*Yas[Aa]#1.#fraction[A]/Xtot
    
#    filename = '../table/'+str(A)+'.txt'
#    filename2 = 'heat'+str(A)+'.dat'
    
#A, Z, Q[MeV], Egamma[MeV], Eelec[MeV], Eneutrino[MeV], tau[s]
#    fchain = pd.read_csv(filename,delim_whitespace=True,header=None)
#    print 'length of the chain',A,len(fchain)

                tmp = []
    
#    for i in range(0,NN):
#        tmp.append(1.0/fchain[6][i])
#    lambdas = np.array(tmp)
####determine the thermalization time in units of day for each element
                tes = np.zeros(NN)
                total_numb = np.zeros(NN)
                
                tas = np.zeros(NN)
                total_numb_a = np.zeros(NN)
    
                for i in range(0,NN):
                    Z = fchain[1][i]
                    Ealpha = Q_alphas[kk][i]
                    Egamma = Q_gammas[kk][i]
                    Eele = Q_betas[kk][i]
       
                    if(Eele>0.):
                        tes[i] = th.calc_thermalization_time(Eele,Mej,vej,Aave,alpha_max,alpha_min,n)
                    if(Ealpha>0.):
                        tas[i] = th.calc_thermalization_time_alpha(Ealpha,Mej,vej,Aave,alpha_max,alpha_min,n)
    #ts = []
                tmp_numb = np.zeros(NN)
                tmp_numb_a = np.zeros(NN)
                number10s = []
                number11s = []
                number12s = []
                number13s = []
                number14s = []
                number15s = []

                heats = []
                gammas = []
                elects = []
                alphas = []

                lambda_sort = np.zeros((NN,NN))
                coeffs = np.ones(NN)
                xs = np.zeros((NN,NN))

   

                for i in range(0,NN):
        
                    tmp0 = lambdas[kk][:i+1]
                    tmp = np.sort(tmp0)
    
                    lambda_sort[i][:i+1] = tmp[::-1]

                for j in range(1,NN):
                    for i in range(0,j):
                        coeffs[j] *= lambda_sort[j-1][i]
   # print j,coeffs[j]

                for j in range(1,NN):
                    for i in range(0,j):
                        xs[j][i] = (lambda_sort[j][j-1-i]-lambda_sort[j][j])
    
    
                for k in range(0,len(ts)):
   # while t<t_final:
                    t = ts[k]
        
                    for i in range(0,NN):
                        coeff = coeffs[i]*np.power(t,i)*branch_chains[kk][i]
        
                        if(i==10):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_10(xs[10][0]*t,xs[10][1]*t,xs[10][2]*t,xs[10][3]*t,xs[10][4]*t,xs[10][5]*t,xs[10][6]*t,xs[10][7]*t,xs[10][8]*t,xs[10][9]*t)
                        elif(i==9):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_9(xs[9][0]*t,xs[9][1]*t,xs[9][2]*t,xs[9][3]*t,xs[9][4]*t,xs[9][5]*t,xs[9][6]*t,xs[9][7]*t,xs[9][8]*t)
                        elif(i==8):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_8(xs[8][0]*t,xs[8][1]*t,xs[8][2]*t,xs[8][3]*t,xs[8][4]*t,xs[8][5]*t,xs[8][6]*t,xs[8][7]*t)
                        elif(i==7):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_7(xs[7][0]*t,xs[7][1]*t,xs[7][2]*t,xs[7][3]*t,xs[7][4]*t,xs[7][5]*t,xs[7][6]*t)
                        elif(i==6):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_6(xs[6][0]*t,xs[6][1]*t,xs[6][2]*t,xs[6][3]*t,xs[6][4]*t,xs[6][5]*t)
                        elif(i==5):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_5(xs[5][0]*t,xs[5][1]*t,xs[5][2]*t,xs[5][3]*t,xs[5][4]*t)          
                        elif(i==4):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_4(xs[4][0]*t,xs[4][1]*t,xs[4][2]*t,xs[4][3]*t)
                        elif(i==3):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_3(xs[3][0]*t,xs[3][1]*t,xs[3][2]*t)
                        elif(i==2):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_2(xs[2][0]*t,xs[2][1]*t)
                        elif(i==1):
                            tmp_numb[i] = coeff*np.exp(-t*lambda_sort[i][i])*bt.calc_M0_1(xs[1][0]*t)
                        elif(i==0):
                            tmp_numb[i] = np.exp(-t*lambda_sort[i][i])
                        else:
                            print 'chain is too long'
         
            
        
    
                    heat = 0.0
                    gam = 0.0
                    ele = 0.0
                    alpha = 0.0
                    ele_th = 0.0
                    gam_th = 0.0
                    alpha_th = 0.0
        
                    for i in range(0,NN):
            
                        Eele = Q_betas[kk][i]
                        Ealpha = Q_alphas[kk][i]
                        Egamma = Q_gammas[kk][i]
            #thermalization for electrons
                        if(t > 0.003*tes[i]):
                            if(Eele > 0.):
                    #tau1 = t/tes[i]
                    #tau0 = 0.05*tau1
                                tau1 = t/tes[i]
                                if(tau1<2.):
                                    tau0 = 0.03*tau1#0.05*tau1                                                                                                                        
                                    root = optimize.newton(th.calc_zero_energy, tau0, args=(tau1,Eele,))
                                    tau0 = root
                                else:
                                    tau0 = 0.4#0.05*tau1  
                    #root = optimize.newton(th.calc_zero_energy, tau0, args=(tau1,Eele,))
                    #tau0 = root
                                delta_t = (tau1-tau0)*tes[i]
            #fth = fchain[6][i]*tes[i]*tes[i]*(np.exp(delta_t/fchain[6][i])-1.0)*np.power(t,-3.)
           # print t,delta_t,delta_t/t,fchain[6][i]*tes[i]*tes[i]*(np.exp(delta_t/fchain[6][i])-1.0)*np.power(t,-3.)
                                t_th = tau0*tes[i]
#                dt_th = delta_t/float(Nth-1)
                                tmp_n_t_th = 1./float(Nth-1)
                                dt_th = np.power(tau1/tau0,tmp_n_t_th)-1.
                                total_numb[i] = 0.0
                                for j in range(0,Nth):
                                    coeff = coeffs[i]*np.power(t_th,i)*branch_chains[kk][i]
                                    tau = t_th/tes[i]
                                    e_delay= th.epsilon_tau(tau,tau1,Eele)#th.calc_e_tau_tau0(tau,tau1)
                        
                                    if(i==10):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_10(xs[10][0]*t_th,xs[10][1]*t_th,xs[10][2]*t_th,xs[10][3]*t_th,xs[10][4]*t_th,xs[10][5]*t_th,xs[10][6]*t_th,xs[10][7]*t_th,xs[10][8]*t_th,xs[10][9]*t_th)
                                    elif(i==9):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_9(xs[9][0]*t_th,xs[9][1]*t_th,xs[9][2]*t_th,xs[9][3]*t_th,xs[9][4]*t_th,xs[9][5]*t_th,xs[9][6]*t_th,xs[9][7]*t_th,xs[9][8]*t_th)
                                    elif(i==8):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_8(xs[8][0]*t_th,xs[8][1]*t_th,xs[8][2]*t_th,xs[8][3]*t_th,xs[8][4]*t_th,xs[8][5]*t_th,xs[8][6]*t_th,xs[8][7]*t_th)            
                                    elif(i==7):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_7(xs[7][0]*t_th,xs[7][1]*t_th,xs[7][2]*t_th,xs[7][3]*t_th,xs[7][4]*t_th,xs[7][5]*t_th,xs[7][6]*t_th)
                                    elif(i==6):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_6(xs[6][0]*t_th,xs[6][1]*t_th,xs[6][2]*t_th,xs[6][3]*t_th,xs[6][4]*t_th,xs[6][5]*t_th)
                                    elif(i==5):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_5(xs[5][0]*t_th,xs[5][1]*t_th,xs[5][2]*t_th,xs[5][3]*t_th,xs[5][4]*t_th)
                                    elif(i==4):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_4(xs[4][0]*t_th,xs[4][1]*t_th,xs[4][2]*t_th,xs[4][3]*t_th)
                                    elif(i==3):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_3(xs[3][0]*t_th,xs[3][1]*t_th,xs[3][2]*t_th)
                                    elif(i==2):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_2(xs[2][0]*t_th,xs[2][1]*t_th)
                                    elif(i==1):
                                        tmp_numb[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_1(xs[1][0]*t_th)
                                    elif(i==0):
                                        tmp_numb[i] = e_delay*np.exp(-t_th*lambda_sort[i][i])                                                                                  
                        
                                    total_numb[i] += tmp_numb[i]*dt_th*t_th
                                                                                  
#                    t_th += dt_th
                                    t_th = (1.+dt_th)*t_th
                                if(Egamma>0.):
                                    Z = fchain[1][i]
                                    kappa_eff = kappa_effs[A][Z]
                                    fth_gamma = th.calc_gamma_deposition(kappa_eff,t,Mej,vej,alpha_min,alpha_max,n) 
              #  print A, Z, kappa_eff, fth_gamma
            #print "tp: ", t/day, A, Z, Egamma,fchain[3][i],Eele,Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A)),Xfraction*MeV*tmp_numb[i]*Eele*lambdas[i]/(mu*float(A))
                                heat += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*(br_betas[kk][i]*Q_total_betas[kk][i])*lambdas[kk][i]/(mu*float(Aa))
                                gam += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_gammas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                                gam_th += fraction_of_chain[kk][i]*fth_gamma*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_gammas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                                ele += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_betas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                                ele_th += fraction_of_chain[kk][i]*np.power(tes[i],2.)*np.power(t,-3.)*Xfraction*MeV*total_numb[i]*br_betas[kk][i]*Q_betas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
          
                        else:
                            if(Egamma>0.):
                                Z = fchain[1][i]
                                kappa_eff = kappa_effs[A][Z]
                                fth_gamma = th.calc_gamma_deposition(kappa_eff,t,Mej,vej,alpha_min,alpha_max,n) 
              #  print A, Z, kappa_eff, fth_gamma
            #print "tp: ", t/day, A, Z, Egamma,fchain[3][i],Eele,Xfraction*MeV*tmp_numb[i]*fchain[3][i]*lambdas[i]/(mu*float(A)),Xfraction*MeV*tmp_numb[i]*Eele*lambdas[i]/(mu*float(A))
                            heat += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*(br_betas[kk][i]*Q_total_betas[kk][i])*lambdas[kk][i]/(mu*float(Aa))
                            gam += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_gammas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                            gam_th += fraction_of_chain[kk][i]*fth_gamma*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_gammas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                            ele += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_betas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                            ele_th += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_betas[kk][i]*Q_betas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
          
                
            #thermalization for alpha particles
                        if(t > 0.003*tas[i]):
                            if(Ealpha > 0.):
                                tau1 = t/tas[i]
                    #tau0 = 0.05*tau1
                    #root = optimize.newton(calc_zero_energy_alpha, tau0, args=(tau1,))
                    #tau0 = root
                                if(tau1<2.):
                                    tau0 = 0.03*tau1#0.05*tau1                                                                                                                        
                                    root = optimize.newton(th.calc_zero_energy_alpha, tau0, args=(tau1,))
                                    tau0 = root
                                else:
                                    tau0 = 0.4#0.05*tau1  
                                delta_t = (tau1-tau0)*tas[i]
                                t_th = tau0*tas[i]
#                dt_th = delta_t/float(Nth-1)
                                tmp_n_t_th = 1./float(Nth-1)
                                dt_th = np.power(tau1/tau0,tmp_n_t_th)-1.
                                total_numb_a[i] = 0.0
                                for j in range(0,Nth):
                                    coeff = coeffs[i]*np.power(t_th,i)*branch_chains[kk][i]
                                    tau = t_th/tas[i]
#                        e_delay= th.epsilon_tau(tau,tau1,Ealpha)
                                    e_delay= th.epsilon_tau_alpha(tau,tau1,Ealpha)
                        #e_delay = calc_e_tau_tau0(tau,tau1)
                    
                                    if(i==10):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_10(xs[10][0]*t_th,xs[10][1]*t_th,xs[10][2]*t_th,xs[10][3]*t_th,xs[10][4]*t_th,xs[10][5]*t_th,xs[10][6]*t_th,xs[10][7]*t_th,xs[10][8]*t_th,xs[10][9]*t_th)
                                    elif(i==9):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_9(xs[9][0]*t_th,xs[9][1]*t_th,xs[9][2]*t_th,xs[9][3]*t_th,xs[9][4]*t_th,xs[9][5]*t_th,xs[9][6]*t_th,xs[9][7]*t_th,xs[9][8]*t_th)
                                    elif(i==8):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_8(xs[8][0]*t_th,xs[8][1]*t_th,xs[8][2]*t_th,xs[8][3]*t_th,xs[8][4]*t_th,xs[8][5]*t_th,xs[8][6]*t_th,xs[8][7]*t_th)            
                                    elif(i==7):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_7(xs[7][0]*t_th,xs[7][1]*t_th,xs[7][2]*t_th,xs[7][3]*t_th,xs[7][4]*t_th,xs[7][5]*t_th,xs[7][6]*t_th)
                                    elif(i==6):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_6(xs[6][0]*t_th,xs[6][1]*t_th,xs[6][2]*t_th,xs[6][3]*t_th,xs[6][4]*t_th,xs[6][5]*t_th)
                                    elif(i==5):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_5(xs[5][0]*t_th,xs[5][1]*t_th,xs[5][2]*t_th,xs[5][3]*t_th,xs[5][4]*t_th)
                                    elif(i==4):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_4(xs[4][0]*t_th,xs[4][1]*t_th,xs[4][2]*t_th,xs[4][3]*t_th)
                                    elif(i==3):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_3(xs[3][0]*t_th,xs[3][1]*t_th,xs[3][2]*t_th)
                                    elif(i==2):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_2(xs[2][0]*t_th,xs[2][1]*t_th)
                                    elif(i==1):
                                        tmp_numb_a[i] = e_delay*coeff*np.exp(-t_th*lambda_sort[i][i])*bt.calc_M0_1(xs[1][0]*t_th)
                                    elif(i==0):
                                        tmp_numb_a[i] = e_delay*np.exp(-t_th*lambda_sort[i][i])                                                                                  
                    
                   
                                    total_numb_a[i] += tmp_numb_a[i]*dt_th*t_th
                                                                                  
                    #t_th += dt_th
                                    t_th = (1.+dt_th)*t_th
            #thermalization for gamma-rays
                                heat += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb_a[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                                alpha += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb_a[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                                alpha_th += fraction_of_chain[kk][i]*np.power(tas[i],2.)*np.power(t,-3.)*Xfraction*MeV*total_numb_a[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                        else:
                            heat += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                            alpha += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
                            alpha_th += fraction_of_chain[kk][i]*Xfraction*MeV*tmp_numb[i]*br_alphas[kk][i]*Q_alphas[kk][i]*lambdas[kk][i]/(mu*float(Aa))
    

   

        
        
                    total_heats[k] += heat
                    total_gammas[k] += gam
                    total_elects[k] += ele
                    total_alphas[k] += alpha
                    total_elect_ths[k] += ele_th
                    total_gamma_ths[k] += gam_th
                    total_alpha_ths[k] += alpha_th
        
                    each_heats[k] += heat
                    each_gammas[k] += gam
                    each_elects[k] += ele
                    each_alphas[k] += alpha
                    each_elect_ths[k] += ele_th
                    each_gamma_ths[k] += gam_th
                    each_alpha_ths[k] += alpha_th
                #print ts[k]/day
#    if(output_each_flag==1):

    data = {'t':np.multiply(ts,1./day),'total':total_heats,'alpha':total_alphas,'gamma':total_gammas, 'electron':total_elects,'alpha_th':total_alpha_ths, 'gamma_th':total_gamma_ths,'electron_th':total_elect_ths}

    print 'end'
    return data

   





def hokan(x,x1,x2,y1,y2):
    return (x2-x)*y1/(x2-x1) + (x-x1)*y2/(x2-x1)
#Katz integral for a given final time
def calc_Katz_integral(Mej,data,t_kf):
    ts = data['t']*day
    integs = []
    integ = ts[0]*ts[0]*(data['gamma_th'][0]+data['electron_th'][0])/(day*day)
    integs.append(integ)
    for i in range(0,len(ts)-1):
        dt = ts[i+1] - ts[i]
        integ += 0.5*(ts[i+1]*data['gamma_th'][i+1]+ts[i+1]*data['electron_th'][i+1] + ts[i]*data['gamma_th'][i]+ts[i]*data['electron_th'][i])*dt/(day*day)
        integs.append(integ)
    k = 0
    while(data['t'][k]<t_kf):
        k +=1
    katz = hokan(t_kf,data['t'][k-1],data['t'][k],integs[k-1],integs[k])
    return katz*Mej
#Katz integral time series
def calc_Katz_integral_timeseries(Mej,data):
    ts = data['t']*day
    integs = []
    integ = ts[0]*ts[0]*(data['gamma_th'][0]+data['electron_th'][0])/(day*day)

    integs.append(integ)

    for i in range(0,len(ts)-1):
        dt = ts[i+1] - ts[i]
        integ += 0.5*(ts[i+1]*data['gamma_th'][i+1]+ts[i+1]*data['electron_th'][i+1] + ts[i]*data['gamma_th'][i]+ts[i]*data['electron_th'][i])*dt/(day*day)
        integs.append(integ*Mej)
    output = {'t':data['t'],'katz':integs}


    return output
