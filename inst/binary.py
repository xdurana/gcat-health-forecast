# Crear datasets binari i categòric
## Fitxers d'entrada

#inst/extdata/gcat/data.csv
#inst/extdata/gcat/categorical.txt -> variables categòriques que convertirà en binari
#inst/extdata/gcat/remove.txt -> variables que no incorporarà al datasets

## Fitxers de sortida

#output/gcat/binary.csv
#output/gcat/categorical.csv

## Script
#python binary.py

# coding: utf-8

# # Build binary dataset

# In[1]:

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import os


# ## Read data and configuration files

# In[17]:

categorical = pd.read_csv('extdata/gcat/categorical.txt', header=None)
remove = pd.read_csv('extdata/gcat/remove.txt', header=None)
data = pd.read_csv('extdata/gcat/data.csv')


# ## Set nulls

# In[18]:

df = data.drop(remove[0].values, axis=1)

for variable in categorical[0].values:
    if variable in df:
        df[df[variable] == 0].variable = None


# In[19]:

df[df.columns[df.columns.str.contains('ETNIA_PARTICIPANTE')]].head()


# ## Write categorical dataset

# In[20]:

df.to_csv('../output/gcat/categorical.csv', index=False)


# ## Transform categorical variables to binary

# In[21]:

for variable in categorical[0].values:
    if variable in df:        
        df[variable] = df[variable].astype('category')
        dummies = pd.get_dummies(df[variable])
        dummies = dummies.rename(columns=lambda x: str('%s_%s' % (variable, int(x))))
        df = pd.concat([df, dummies], axis=1)
        df = df.drop(variable, axis=1)


# In[22]:

df[df.columns[df.columns.str.contains('ETNIA_PARTICIPANTE_1') |
              df.columns.str.contains('ETNIA_PARTICIPANTE_2') |
              df.columns.str.contains('ETNIA_PARTICIPANTE_3') |
              df.columns.str.contains('ETNIA_PARTICIPANTE_4')]].head()


# ## Write binary dataset

# In[23]:

df.to_csv('../output/gcat/binary.csv', index=False)
