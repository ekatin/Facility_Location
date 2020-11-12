#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import random

import gurobipy as gp
from gurobipy import GRB

import matplotlib.pyplot as plt

from itertools import product
from math import sqrt

import numpy as np


# In[16]:


ilceler = [(40.874572, 29.132615),(41.185826, 28.720197),(40.997592, 29.100956),(40.984296, 28.723891),(41.033995, 28.832999),
          (40.998019, 28.851395),(40.981307, 28.873582),(41.118848, 28.804236),(41.041590, 28.908753),(41.074971, 29.020352),
          (41.126921, 29.099183),(41.011481, 28.650702),(41.028456, 28.973512),(41.021081, 28.583458),(41.142084, 28.457303),
          (41.033176, 29.168299),(41.054168, 28.867725),(41.019303, 28.687639),(41.043394, 28.929759),(41.017764, 28.940846),
          (41.084140, 28.888458),(41.010768, 28.874884),(40.993689, 29.037498),(41.072188, 28.964582),(40.890213, 29.183889),
          (41.035418, 28.785618),(40.931593, 29.128047),(40.891644, 29.264486),(40.991552, 29.230045),(41.111961, 29.027793),
          (41.072951, 28.255943),(40.970196, 29.260874),(41.106731, 28.872237),(41.177618, 29.611007),(41.064482, 28.983654),
          (40.843431, 29.300695),(41.029830, 29.099501),(41.020928, 29.018995),(40.996549, 28.909799)]

adaylar = [(40.96953784, 29.05790448),(40.93056500, 29.11498400),(40.95333500, 29.27768800),(40.95424967, 29.08190489),(40.96493637, 29.05599475),
           (40.89036100, 29.25383900),(41.03164100, 29.02594900),(41.02586000, 29.06967600),(41.01149700, 29.16110500),(41.13148800, 29.09446400),
           (41.10483000, 29.07460000),(41.10932500, 29.05177300),(41.04911500, 29.01493300),(41.1806430, 28.74963100),(40.97098383, 28.73813689),
           (40.97061930, 28.72102439),(40.97030900, 28.79192800),(40.95784600, 28.81205500),(41.10756607, 28.67121041),(41.09286380,28.64975542),
           (41.03862900, 28.98667300),(41.05614300, 28.94912100),(41.00101300, 28.94441500),(41.06717600,28.95821200),(40.98190985,28.76146711),
           (41.00007900, 28.76459900),(41.09560300, 28.89952700),(41.04344500, 28.99427900),(41.02187100, 28.92142800)]



# In[4]:


max_facilities=10

Kapsama_mesafesi = 25

P=[0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0]
# 1. Aşama sonucunla bulunan noktaların koordinatları
P_koordinatlari=[(40.998019, 28.851395),(41.118848, 28.804236),(41.054168, 28.867725),(41.019303, 28.687639),
                 (41.072188, 28.964582),(41.035418, 28.785618),(40.891644, 29.264486),(40.991552, 29.230045),
                 (41.02983, 29.099501),(41.020928, 29.018995)]


# In[7]:


R = 6371

def pol2cart(a, i):
    a, i = np.radians(a), np.radians(i)
    return R*np.cos(a) *np.cos(i),           R*np.cos(a) *np.sin(i),           R*np.sin(a)        
      
def haversine_dist(a,i): 
    point1_cart = np.array(pol2cart(*a))
    point2_cart = np.array(pol2cart(*i))
    euc_dist = np.linalg.norm(point1_cart-point2_cart)
    sin_theta_2 = euc_dist / (R * 2)     
    theta_2 = np.arcsin(sin_theta_2)  
    theta = 2*theta_2
    dist = R*theta
    return dist
    
for a in adaylar:
    for i in ilceler:
         haversine_dist(a,i)
    
            
dist_ai=[]

#print(len(dist_ai)) 

for a in adaylar:
    for i in ilceler:
           dist_ai.append(haversine_dist(a,i))

#print(len(dist_ai))    
#np.mean(dist_ai)


# In[8]:


num_ilceler = len(ilceler)
num_adaylar= len(adaylar)

cartesian_prod1 = list(product(range(num_adaylar), range(num_ilceler)))
mesafeler1={(a,i): haversine_dist(adaylar[a], ilceler[i]) for a, i in cartesian_prod1}
                     
eslesmeler1 = {(a,i): haversine_dist(adaylar[a], ilceler[i])
            for a in range(num_adaylar)
            for i in range(num_ilceler) 
            if  haversine_dist(adaylar[a],ilceler[i]) < Kapsama_mesafesi}
print("Geçerli eşleşme sayısı: {0}".format(len(eslesmeler1.keys())))
#print(mesafeler1)
#print(eslesmeler1)


# In[9]:


m = gp.Model("2.asama")

# Karar Değişkenleri : j tesisi açılırsa 1 , dd : 0
select = m.addVars(range(num_adaylar), vtype=GRB.BINARY, name='oj')
# Karar Değişkenleri: k talep noktası, j tesisinden hizmet alırsa:1 , dd:0
assign = m.addVars(eslesmeler1.keys(), vtype=GRB.BINARY, name='Xjk')

# Amaç Fonskyonu
# Min= Djk*Xjk
#m.setObjective(assign.prod(eslesmeler1), GRB.MINIMIZE)
obj = gp.quicksum( eslesmeler1[a,i] *assign[a,i] for a,i in eslesmeler1.keys())
m.setObjective(obj, GRB.MINIMIZE)

# Kısıt 1 
m.addConstr(select.sum() >= max_facilities, name="Top_oj>=F")

# Kısıt 2
m.addConstrs((assign[a,i] == select[a]
             for a, i in eslesmeler1.keys()),
             name="Xj,k=oj")

# Kısıt 3

m.addConstrs((assign.sum(adaylar,'*') <= P[i]
             for i in range(num_ilceler)),
            name="Topj_Xj,k<=Pk")

# Optimum çözüm
m.optimize()


# In[10]:


for v in m.getVars():
        print('%s %g' % (v.varName, v.x))
       # print('Obj: %g' % m.objVal)


# In[17]:


sonuclar=[(40.953335, 29.277688),(40.890361, 29.253839),(40.97098383, 28.73813689),(40.9706193, 28.72102439),
          (40.970309, 28.791928),(41.0928638, 28.64975542),(41.056143, 28.949121),(41.067176, 28.958212),
         (41.000079, 28.764599),(41.095603, 28.899527)]

heliportlar=[(41.0726,29.01172),(41.05674,28.98806),(41.10971,29.01834),(40.99217,29.10645),(41.02161,28.65726),(40.98333,29.05435),
             (40.97118,29.09865),(41.02563,28.977),(41.06796,28.67552),(41.10259,28.98675),(40.9725,28.87092),(41.02701,29.10626),
             (40.9119,29.17888),(41.10676,28.80304),(41.04202,28.83356),(41.00074,28.62727)]

# In[18]:


plt.figure(figsize=(6,6), dpi=110)
plt.scatter(*zip(*sonuclar), c='Blue', s=100)
plt.scatter(*zip(*ilceler), c='Red', s=20)
assignments = [p for p in eslesmeler1 if assign[p].x > 0.5]
for p in assignments:
    pts = [adaylar[p[0]],ilceler[p[1]]]
    #plt.plot(*zip(*pts), c='Black', linewidth=0.1)


# In[21]:


import folium

#map = folium.Map(location=[41,28], zoom_start = 10)

markers_dict = {'Sultanbeyli Semt Parki':[40.953335, 29.277688],
                'Gozdagi Korusu':[40.890361, 29.253839],
                'Avcılar Sahil Parkı-2':[40.97098383, 28.73813689],
                'Avcilar Sahil Parki1':[40.9706193, 28.72102439],
                'Florya Sahil Parki':[40.970309, 28.791928],
                'Toki Hosdere Hayat Parki':[41.0928638, 28.64975542],
                'Sütlüce–2 Parki':[41.056143, 28.949121],
                'Sadabad Mesire –1':[41.067176, 28.958212],
                'Küçükçekmece Sahil Parki':[41.000079, 28.764599],
                'Gazi Parki':[41.095603, 28.899527]}



map_cities= folium.Map(location=[41,28], zoom_start = 8)

for i in markers_dict.items():
    folium.Marker(location=i[1], popup=i[0]).add_to(map_cities)
    print(i)
    
map_cities

