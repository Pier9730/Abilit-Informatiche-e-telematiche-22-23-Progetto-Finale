import astropy.io
from astropy.io import fits
import numpy as np
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
import os

tstart=time.time()


test=1                                                                                                                                                  #1, 2 or 3 
arr_p=[1, 3, 5]                                                                                                                                         #1, 3, or 5 (multipole)
Nsim=10000                                                                                                                                                #number of simulation (FITS files)
Nval=200                                                                                                                                                #number of values for each simulation
debug=0                                                                                                                                                 #0=no debug; 1=first level; 2=print diagonals; 3=print all matrix el; 4=all corr el; 5=print ave_list

directory= "Test%i_Nval_%i_Nsim_%i_corr" %(test, Nval, Nsim)
home_dir="/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/Immagini"
path=os.path.join(home_dir, directory)
os.mkdir(path)


if test==1:
    sig_l=[0, 0.02, 0, 0.02, 0, 0.02]
    h_l=[0, 25, 0, 50, 0, 75]
elif test==2:
    sig_l=[0, 0.02, 0, 0.01, 0, 0.005]
    h_l=[0, 50, 0, 50, 0, 50]
else:
    sig_l=[0, 0.02, 0, 0.01, 0, 0.005]
    h_l=[0, 5, 0, 5, 0, 5]
    


data_list=[] 
for k in range(Nsim):
    hdul=fits.open('/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/data/MockMeasures_2PCF_Test%i/MockMeasures_2PCF_Correlation_MULTIPOLES_Test%i_%i.fits' %(test, test, k+1))
    #hdul.info()
    data_list.append(hdul[1].data)
    hdul.close()
    
pos_list=[]    
for i in range(Nval):                                                                                                                                   #list of positions on the grid
    pos_list.append(data_list[0][i][0])

pole_list=[]                                                                                                                                            #list of multipole values
for p in range(1,6):
    for i in range(Nval):
        for k in range(Nsim):
            pole_list.append(data_list[k][i][p])

print("pole_list lenght:", len(pole_list))    



if debug==1:
    print("data_list lenght:", len(data_list))
    print(data_list[0][0][0])                                                                                                                           #first index: fits file
                                                                                                                                                        #second index: which of 200 6-uple
                                                                                                                                                        #third index: which element of the 6-uple; 2 and 4 useless (even multipoles)
if debug==1:
                                               
    temp=[] 
    temp2=[]                                                   
    for i in range(Nsim):
        temp.append(data_list[i][0][arr_p[0]])
        temp2.append((data_list[i][0][arr_p[0]])**2)
    sigma_2=np.average(temp2)-np.average(temp)**2
    print("Var of first monopole:", sigma_2)

#COMPUTE NUMERICAL AND THEORETICAL COVARIANCE MATRIX

def cov_many(x,y,sig1, h1, sig2, h2):                                                                                                                   #covariance with correlation between multipoles
    return np.exp(-(x-y)**2/(h1**2+h2**2))*np.sqrt(2*h1*h2/(h1**2+h2**2))*sig1*sig2




cov_num=np.zeros((3*Nval,3*Nval),dtype=float)
cov_th=np.zeros((3*Nval,3*Nval),dtype=float)
res=np.zeros((3*Nval,3*Nval),dtype=float) 
corr_th=np.zeros((3*Nval,3*Nval),dtype=float)
 
ave_list=[]
for p in arr_p: 
    for i in range(Nval):
        vallist=[]
        for j in range(Nsim):
            vallist.append(data_list[j][i][p])
        ave_list.append(np.average(vallist)) 
                                                                                                                                                        #appends i-th average on the list
for p1 in arr_p:
    for p2 in arr_p:
        print("p1 index:", arr_p.index(p1))
        print("p2 index:", arr_p.index(p2))
        
        p2_shift=Nval*arr_p.index(p2)
        p1_shift=Nval*arr_p.index(p1)
        
        for i in range(Nval):                                                                
            if i%10==0:
                print(i)
            for j in range(Nval):
                cov_th[p2_shift+i, p1_shift+j]=cov_many(pos_list[i],pos_list[j], sig_l[p1], h_l[p1], sig_l[p2], h_l[p2])                                #theoretical covariance
                #if debug==3:
                #    print("r_i and r_j:",data_list[0][i][0], " || ", data_list[0][j][0], " ||| cov_th[i,j]=",cov_th[i,j] )
                summ=0
                for k in range(Nsim):
                    #summ+=(data_list[k][i][p2]-ave_list[p2_shift+i])*(data_list[k][j][p1]-ave_list[p1_shift+j])                                         #numerical covariance 
                    summ+=(pole_list[(p2-1)*Nval*Nsim+i*Nsim+k]-ave_list[p2_shift+i])*(pole_list[(p1-1)*Nval*Nsim+j*Nsim+k]-ave_list[p1_shift+j])                             #DA METTERE A POSTO                        
                cov_num[p2_shift+i, p1_shift+j]=summ/(Nsim-1) 
        """
        print("cov_num:", cov_num)
        print("cov_th:", cov_th)
        
        plt.figure(1)
        plt.title('Covariance num ')
        plt.imshow(cov_num)
        #plt.axis('off')
        plt.show()

        plt.figure(2)
        plt.title('Covariance th')
        plt.imshow(cov_th)
        #plt.axis('off')
        plt.show()                                                                                     
        """

        if debug==1: 
    
            print("Covariance:", type(cov_num))  
        if debug==2:
            diag=[]
            diag_th=[]
            for i in range(Nval):
                for j in range(Nval):
                    if j==i:
                        diag.append(cov_num[i,j]) 
                        diag_th.append(cov_th[i,j]) 
            print("Diagonal:", diag) 
            
            print("Theoretical diagonal:", diag_th)                                                                              
      
    
    


        
#COMPUTE CORRELATION AND RESIDUAL MATRIX

  
        if debug==2:
            diag_corr=[]                                                   
        for i in range(Nval):
            for j in range(Nval):
                corr_th[p2_shift+i,p1_shift+j]=cov_th[p2_shift+i,p1_shift+j]/np.sqrt(cov_th[p2_shift+i,p1_shift+i]*cov_th[p2_shift+j,p1_shift+j])
                corr2=corr_th[p2_shift+i,p1_shift+j]
                if debug==4:
                    print("i:",i, " || j :",j, " ||| corr_th[i,j]=",cov_th[i,j] )
                if debug==2 and i==j:
                    diag_corr.append(corr_th[i,j])            
                res[p2_shift+i,p1_shift+j]=(cov_th[p2_shift+i,p1_shift+j]-cov_num[p2_shift+i,p1_shift+j])*np.sqrt((Nsim-1)/((1+corr2**2)*cov_th[p2_shift+i,p1_shift+i]*cov_th[p2_shift+j,p1_shift+j]))

            
        if debug==2:
            print("Correlation diagonal:", diag_corr)
        if debug==5:
            print("Cov theory lenght:", len(cov_th))                
            print("covariance lenght:", len(cov_num))
            print("Corr theory lenght:", len(corr_th))                
            print("residual lenght:", len(res))

cov_num=cov_num/(Nsim-1)                                                                                                                                #divide every element of cov_num by Nsim -1        
#COMPUTE STANDARD DEVIATION OF RESIDUAL
res_arr=res.reshape((3*Nval)**2)
res_dev=np.std(res_arr)
print("Standard deviation:", res_dev)                                                        
    

if res_dev<1.1:
    print("!!!!!!!!!!!!!!\n      =)      \n!!!!!!!!!!!!!!")
else:
    print("!!!!!!!!!!!!!!\n      =(      \n!!!!!!!!!!!!!!")

tstop=time.time()
deltat=tstop-tstart
print("Elapsed time:", deltat)           



fig=plt.figure(1)
plt.title('Covariance num')
plt.imshow(cov_num)
plt.show()
fig.savefig("C:/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/Immagini/%s/cov_num.jpg" %(directory))

fig=plt.figure(2)
plt.title('Covariance th')
plt.imshow(cov_th)
plt.show()
fig.savefig("C:/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/Immagini/%s/cov_th.jpg" %(directory))

fig=plt.figure(3)
plt.title('Residual')
plt.imshow(res)
plt.show()
fig.savefig("C:/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/Immagini/%s/res.jpg" %(directory))


fig=plt.figure(4)
plt.hist(res_arr, 50, density=True)
plt.show()
fig.savefig("C:/Users/Pierfrancesco/Abilità_Informatiche/Progetto_Finale/Immagini/%s/hist.jpg" %(directory))






