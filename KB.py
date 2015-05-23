#!/usr/bin/python
import numpy as np
import copy as copy
import math
import matplotlib.pyplot as plt
for file in ['Original_FF','minus_6','minus_4','16','12','8']:
  
 if file == 'Original_FF':
    V=65.1065
    #break                      #average volume of simulatin box
 elif file == 'minus_6':
    V=53.2044
    #break
 elif file=='minus_4':
    V=54.3159
    #break
 elif file=='8':
    V=56.7635
    #break
 elif file=='12':
    V=57.5747
    #break
 elif file=='16':
    V=58.5029
    #break
 inp=open('rdf_%s_25.xvg' % file,'r')

 in_name='RDF'

 R=[]
 g=[]

 for line in inp:
    line=line.split()

    R.append(float(line[0]))


    g.append(float(line[1]))
 print 'R'
 print len(R)
 
##############################delta Nsw(r)############################################

 g_m=[]                         #the function of [gsw(r)-1]*r^2

 delta_N=[]                     # the integration of [gsw(r)-1]*r^#

#d_delta_N=[]                   #deviation of delta_N

 N=1056                         #number of water molecule in box
 
 Avg_density = N/V              #molecular density 

 f1=4*math.pi*Avg_density       #constate 
 f3=4*math.pi

# calculating f1*[gsw(r)-1]*r^2, named as g_m

 i=0
 while i<len(R):
        tem1=f1*(g[i]-1)*((R[i])**2)
        g_m.append(tem1)
        i=i+1


# integral of f1*[gsw(r)-1]*r^2, named as delta_N

 j=1
 sum=0
 while j<len(R):
        sum=sum+(g_m[j-1]+g_m[j])*(R[j]-R[j-1])*0.5
        delta_N.append(sum)
        j=j+1

 print 'delta_N'
 print len(delta_N)
 out=open('delta_N_%s.xvg' % file,'w')
 k=0
 while k<len(delta_N):
      out.write(str(R[k])+' '+str(delta_N[k])+'\n')
      k=k+1
 out.close()


##############################corrected rdf########################################

 g_c=[]                        #correction of radail distribution function  

 f2=(4/3)*math.pi              #constant

 k=0
# calculation of corrected radial distribution function, named as g_correct 
 while k<len(R)-1:
        Volume_rate =1-(f2*(R[k]**3)/43.2952)
        tem3=(g_m[k])*N*Volume_rate/((N*Volume_rate)-delta_N[k])
        g_c.append(tem3)
        k=k+1

 print 'g_c' 
 print len(g_c)


################################RKBI calculation#######################################

 G=[]                           # Kirkwood buff integral

 delta_G=[]                     # integrating function [g_c-1]*r^2                

 f3=4*math.pi                   #constant 

#calculation of integrating function [g_c-1]*r^2, named as delta_G
 u=0
 while u<len(g_c):
        tem4=f3*(g_c[u]-1)*((R[u])**2)
        delta_G.append(tem4)
        u=u+1
 print 'delta_G'
 print len(delta_G)
#calculation of RKBI, named as G

 h=1
 sum=0
 while h<len(delta_G):
        sum=sum+(delta_G[h-1]+delta_G[h])*(R[h]-R[h-1])*0.5
        G.append(sum)
        h=h+1

 print 'G' 
 print len(G)


##############################save and plot the corrected RDF#######################

 R_1=[]                        # since the dimension of R and g_c isn't same, let's make a new R_1

 m=0
 while m<len(G):
      R_1.append(R[m])
      m=m+1



 plt.plot(R_1,G,'r-')
 plt.show()




 R_2=[]
 m=0
 while m<len(delta_N):
      R_2.append(R[m])
      m=m+1
 plt.plot(R_2,delta_N,'r-')
 plt.show()









