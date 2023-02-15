import astropy.io
from astropy.io import fits
import numpy as np
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
import os




tstart=time.time()


test=2                                                                                      #1, 2 or 3 
p=1                                                                                         #1, 3, or 5 (multipole)
Nsim=100                                                                                    #number of simulation (FITS files)
Nval=200                                                                                    #number of values for each simulation
debug=0                                                                                     #0=no debug; 1=first level; 2=print diagonals; 3=print all matrix el; 4=all corr el

directory= "Test%i_Nval_%i_Nsim_%i_pole%i" %(test, Nval, Nsim, p)
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
    #for i in range(Nval):
    #     pole_list.append(data_list[k][i][p]) 
    
pole_list=[]                                                                                #list of multipole values
for i in range(Nval):
    for k in range(Nsim):
        pole_list.append(data_list[k][i][p])
        
print(len(pole_list))


pos_list=[]                                                                                 #list of positions on the grid
for i in range(Nval):
    pos_list.append(data_list[0][i][0])
    
tstop1=time.time()
print("load time=", tstop1-tstart)




if debug==1:
    print("data_list lenght:", len(data_list))
    print(data_list[1][0][0])                                                               #first index: fits file
                                                                                            #second index: which of 200 6-uple
                                                                                            #third index: which element of the 6-uple; 2 and 4 useless (even multipoles)
if debug==1:
                                               
    temp=[] 
    temp2=[]                                                   
    for i in range(Nsim):
        temp.append(data_list[i][0][p])
        temp2.append((data_list[i][0][p ])**2)
    sigma_2=np.average(temp2)-np.average(temp)**2
    print("Var of first monopole:", sigma_2)

#COMPUTE NUMERICAL AND THEORETICAL COVARIANCE MATRIX

def cov_sing(x,y,sig, h):                                                                   #covariance without correlation between multipoles
    return sig**2*np.exp(-((x-y)**2)/(2*h**2))



cov_num=np.zeros((Nval,Nval),dtype=float)
cov_th=np.zeros((Nval,Nval),dtype=float)

ave_list=[]
for i in range(Nval):
    vallist=[]
    for j in range(Nsim):
        vallist.append(data_list[j][i][p])
    ave_list.append(np.average(vallist))                                                    #appends i-th average on the list

#print("ave_list lenght:", len(ave_list))
#print("p0 average for first distance:", ave_list[0])


 


for i in range(Nval):
    if i%10==0:
        print(i)
    for j in range(Nval):
        cov_th[i,j]=cov_sing(pos_list[i],pos_list[j],sig_l[p], h_l[p])        #theoretical covariance
        if debug==3:
            print("r_i and r_j:",data_list[0][i][0], " || ", data_list[0][j][0], " ||| cov_th[i,j]=",cov_th[i,j] )

        
        summ=0
        for k in range(Nsim):
            #summ+=(pole_list[k*Nval+i]-ave_list[i])*(pole_list[k*Nval+j]-ave_list[j])         #numerical covariance  
            summ+=(pole_list[i*Nsim+k]-ave_list[i])*(pole_list[j*Nsim+k]-ave_list[j])                                 

        cov_num[i,j]=summ                                                                                       
cov_num=cov_num/(Nsim-1)            

if debug==1: 
    print("Cov theory lenght:", len(cov_th))                
    print("covariance lenght:", len(cov_num))
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

res=np.zeros((Nval,Nval),dtype=float) 
corr_th=np.zeros((Nval,Nval),dtype=float)   
if debug==2:
    diag_corr=[]                                                   
for i in range(Nval):
    for j in range(Nval):
        corr_th[i,j]=cov_th[i,j]/np.sqrt(cov_th[i,i]*cov_th[j,j])
        if debug==4:
            print("i:",i, " || j :",j, " ||| corr_th[i,j]=",cov_th[i,j] )
        if debug==2 and i==j:
            diag_corr.append(corr_th[i,j])            
        res[i,j]=(cov_th[i,j]-cov_num[i,j])*np.sqrt((Nsim-1)/((1+corr_th[i,j]**2)*cov_th[i,i]*cov_th[j,j]))

if debug==2:
    print("Correlation diagonal:", diag_corr)
    
#COMPUTE STANDARD DEVIATION OF RESIDUAL
res_arr=res.reshape(Nval**2)
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





