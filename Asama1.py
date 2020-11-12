#!/usr/bin/env python
# coding: utf-8

# In[19]:


from itertools import product
from math import sqrt

import gurobipy as gp
from gurobipy import GRB

import numpy as np


# In[20]:


talepler=[(40.874572, 29.132615),(41.185826, 28.720197),(40.997592, 29.100956),(40.984296, 28.723891),(41.033995, 28.832999),
          (40.998019, 28.851395),(40.981307, 28.873582),(41.118848, 28.804236),(41.041590, 28.908753),(41.074971, 29.020352),
          (41.126921, 29.099183),(41.011481, 28.650702),(41.028456, 28.973512),(41.021081, 28.583458),(41.142084, 28.457303),
          (41.033176, 29.168299),(41.054168, 28.867725),(41.019303, 28.687639),(41.043394, 28.929759),(41.017764, 28.940846),
          (41.084140, 28.888458),(41.010768, 28.874884),(40.993689, 29.037498),(41.072188, 28.964582),(40.890213, 29.183889),
          (41.035418, 28.785618),(40.931593, 29.128047),(40.891644, 29.264486),(40.991552, 29.230045),(41.111961, 29.027793),
          (41.072951, 28.255943),(40.970196, 29.260874),(41.106731, 28.872237),(41.177618, 29.611007),(41.064482, 28.983654),
          (40.843431, 29.300695),(41.029830, 29.099501),(41.020928, 29.018995),(40.996549, 28.909799)]

tesisler=[(40.874572, 29.132615),(41.185826, 28.720197),(40.997592, 29.100956),(40.984296, 28.723891),(41.033995, 28.832999),
          (40.998019, 28.851395),(40.981307, 28.873582),(41.118848, 28.804236),(41.041590, 28.908753),(41.074971, 29.020352),
          (41.126921, 29.099183),(41.011481, 28.650702),(41.028456, 28.973512),(41.021081, 28.583458),(41.142084, 28.457303),
          (41.033176, 29.168299),(41.054168, 28.867725),(41.019303, 28.687639),(41.043394, 28.929759),(41.017764, 28.940846),
          (41.084140, 28.888458),(41.010768, 28.874884),(40.993689, 29.037498),(41.072188, 28.964582),(40.890213, 29.183889),
          (41.035418, 28.785618),(40.931593, 29.128047),(40.891644, 29.264486),(40.991552, 29.230045),(41.111961, 29.027793),
          (41.072951, 28.255943),(40.970196, 29.260874),(41.106731, 28.872237),(41.177618, 29.611007),(41.064482, 28.983654),
          (40.843431, 29.300695),(41.029830, 29.099501),(41.020928, 29.018995),(40.996549, 28.909799)]

nufus= [(15.238),(282.488),(425.094),(448.882 ),(745.125),(611.059),(229.239 ),(460.259),(274.735),(182.649),(248.260),
       (352.412),(233.323 ),(254.103),(73.718 ),(264.508 ),(450.344 ),(954.579 ),(400.513 ),(443.090 ),(491.962 ),(289.441),
       (482.713 ),(448.025),(470.676 ),(792.821 ),(513.316),(711.894 ),(436.733),(347.214),(193.680),(336.021 ),
       (534.565 ),(37.692 ),(279.817),(267.400),(710.280 ),(531.825 ),(293.574 )]
      


# In[5]:


Toplam_tesis=15519.267000000003/1551.9267000000004
#print(Toplam_tesis)
TN=sum(nufus)
#print(TN)
C=15519.267000000003/10
#print(C)
num_talepler = len(talepler)
num_tesisler = len(tesisler)
cartesian_prod = list(product(range(num_talepler), range(num_tesisler)))
num_nufus=len(nufus)


# In[25]:


import numpy as np
R = 6371

def pol2cart(i, j):
    i, j = np.radians(i), np.radians(j)
    return R*np.cos(i) *np.cos(j),           R*np.cos(i) *np.sin(j),           R*np.sin(i)        
      
def haversine_dist(i,j): 
    point1_cart = np.array(pol2cart(*i))
    point2_cart = np.array(pol2cart(*j))
    euc_dist = np.linalg.norm(point1_cart-point2_cart)
    sin_theta_2 = euc_dist / (R * 2)     
    theta_2 = np.arcsin(sin_theta_2)  
    theta = 2*theta_2
    dist = R*theta
    return dist
    
for i in talepler:
    for j in tesisler:
         haversine_dist(i,j)
    
            
dist_ij=[]

#print(len(dist_ij)) 

for i in talepler:
    for j in tesisler:
           dist_ij.append(haversine_dist(i,j))

#print(len(dist_ij))      
#print(dist_ij)  

uzaklık= {(i,j):haversine_dist(talepler[i], tesisler[j]) for i, j in cartesian_prod} #dist(i,j)

#ilçeler arası ortalama uzaklık 
np.mean(dist_ij)
print(uzaklık)


# In[7]:


m = gp.Model('1.Asama')
sayisi = m.addVars ( range(num_talepler),vtype = GRB.INTEGER,  name = "fi_sayisi") 
assign = m.addVars(cartesian_prod,vtype=GRB.CONTINUOUS , name='Xij')
m.setObjective(assign.prod(uzaklık), GRB.MINIMIZE)
# toplam fi = F kısıtı
m.addConstr(sayisi.sum() ==Toplam_tesis , name="toplamtesiskısıt")
# her bölgenin talebini karşılayan kısıt toplamXij=Ri
m.addConstrs((gp.quicksum(assign[(c,f)] for f in range(num_tesisler)) == nufus[c] for c in range(num_talepler)), name='karsılanan_talep')
#karşılanan talep TopXij<= fiC
m.addConstrs((gp.quicksum(assign[(c,f)] for c in range(num_talepler)) <= sayisi[f]*C for f in range(num_tesisler)), name='kapasite_talep')

m.optimize()


# In[8]:


for sayisi in m.getVars():
        print('%s %g' % (sayisi.varName, sayisi.x))
        #print('Obj: %g' % m.objVal)           


# In[ ]:


Secilen_noktalar=[(40.998019, 28.851395),(41.118848, 28.804236),(41.054168, 28.867725),(41.019303, 28.687639)
                 (41.072188, 28.964582),(41.035418, 28.785618),(40.891644, 29.264486),(40.991552, 29.230045)
                 (41.02983, 29.099501),(41.020928, 29.018995)]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




