# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 15:57:33 2018

@author: PierreLoup
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def plot_extent(filename, month):
    
    data = pd.read_csv(filename, sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)
    
    #Discard missing entries
    data = data.drop(data[data[" extent"] < -9000].index)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Year")
    ax.set_ylabel("Ice extent")

    data.iloc[:,3].plot(label="Ice extent in {}".format(month))
    plt.legend(loc='best')

    plt.savefig('extent_{}.pdf'.format(month))
    plt.close()


def plot_area(filename, month):
    
    data = pd.read_csv(filename, sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)
    
    #Discard missing entries
    data = data.drop(data[data[" extent"] < -9000].index)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Year")
    ax.set_ylabel("Ice area")

    data.iloc[:,4].plot(label="Ice area in {}".format(month))
    plt.legend(loc='best')

    plt.savefig('area_{}.pdf'.format(month))
    plt.close()

def plot_extent_variations(filename, month):
    
    data = pd.read_csv(filename, sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)
    
    #Discard missing entries
    data = data.drop(data[data[" extent"] < -9000].index)
    
    #Difference in ice extent between two consecutive years
    diff = data.iloc[:,3].diff()
    x = np.linspace(1, data.shape[0], data.shape[0])
    plt.plot(x, diff, label="Variation")
    diff.cumsum().plot(label="Cum. Variation")
    plt.legend(loc='best')

    plt.savefig('extent_variation_{}.pdf'.format(month))
    plt.close()
    
def plot_area_variations(filename, month):
    
    data = pd.read_csv(filename, sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)
    
    #Discard missing entries
    data = data.drop(data[data[" extent"] < -9000].index)
    
    #Difference in ice extent between two consecutive years
    diff = data.iloc[:,4].diff()
    x = np.linspace(1, data.shape[0], data.shape[0])
    plt.plot(x, diff, label="Variation")
    diff.cumsum().plot(label="Cum. Variation")
    plt.legend(loc='best')

    plt.savefig('area_variation_{}.pdf'.format(month))
    plt.close()

######################### January measurements ######################

dataJanuary = pd.read_csv('N_01_extent_v3.0.csv', sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)

#Discard missing entries
dataJanuary = dataJanuary.drop(dataJanuary[dataJanuary[" extent"] < -9000].index)

plot_extent("N_01_extent_v3.0.csv", "January")
plot_area("N_01_extent_v3.0.csv", "January")

plot_extent_variations("N_01_extent_v3.0.csv", "January")
plot_area_variations("N_01_extent_v3.0.csv", "January")

print("January measurements:")

#Mean ice extent
print("Mean ice extent (millions of squared km):")
mean = dataJanuary.iloc[:,3].mean()
print(mean)
#Mean ice area
print("Mean ice area (millions of squared km):")
mean = dataJanuary.iloc[:,4].mean()
print(mean)

#Mean extent variation
print("Mean ice extent variation (millions of squared km):")
diff = dataJanuary.iloc[:,3].diff()
print(np.mean(diff))
#That's 20 times the size of Belgium, each year
#Mean extent variation
print("Mean ice area variation (millions of squared km):")
diff = dataJanuary.iloc[:,4].diff()
print(np.mean(diff))


print("")
print("")

######################### June measurements ######################

dataJune = pd.read_csv('N_06_extent_v3.0.csv', sep=",", usecols=[0,2,3,4,5], nrows=39, header=0)

#Discard missing entries
dataJune = dataJune.drop(dataJune[dataJune[" extent"] < -9000].index)

plot_extent("N_06_extent_v3.0.csv", "June")
plot_area("N_06_extent_v3.0.csv", "June")

plot_extent_variations("N_06_extent_v3.0.csv", "June")
plot_area_variations("N_06_extent_v3.0.csv", "June")

print("June measurements:")

#Mean ice extent
print("Mean ice extent (millions of squared km):")
mean = dataJune.iloc[:,3].mean()
print(mean)
#Mean ice area
print("Mean ice area (millions of squared km):")
mean = dataJune.iloc[:,4].mean()
print(mean)

#Mean extent variation
print("Mean ice extent variation (millions of squared km):")
diff = dataJune.iloc[:,3].diff()
print(np.mean(diff))

#Mean extent variation
print("Mean ice area variation (millions of squared km):")
diff = dataJune.iloc[:,4].diff()
print(np.mean(diff))