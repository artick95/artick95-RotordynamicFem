import math
import cmath
import matplotlib.pyplot
import random
import turtle
import time


print('how many beam elements we have?')
n_elements=float(input())
l=[1]*n_elements

u=[1]*n_elements*2
teta=[1]*n_elements*2

E=2e09
print('insert the moment of inertia of the area cross section?')
I=float(input())


for i in range(1,n_elements):
    print('insert the '+ str(i)+' elements lenght \n')
    l[i]=float(input())

