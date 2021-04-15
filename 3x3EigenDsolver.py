import math


print('insert d11[m]')
d11=float(input())
print('insert d12[m]')
d12=float(input())
print('insert d13[m]')
d13=float(input())
print('insert d21[m]')
d21=float(input())
print('insert d22[m]')
d22=float(input())
print('insert d23[m]')
d23=float(input())
print('insert d31[m]')
d31=float(input())
print('insert d32[m]')
d32=float(input())
print('insert d33[m]')
d33=float(input())


print('-w^6'+'+w^4*'+str((d11+d22+d33))+'+z^2*'+str(d11*d22-d12*d21+d11*d33-d13*d31+d22*d33-d23*d32)+'+'+str(-d13*d21*d32-d12*d23*d31-d11*d22*d33+d13*d22*d31+d12*d21*d33+d11*d23*d32))

sol_eigen=[]

for eigenvalue in range(0,10**5,0.5):
    eqn=eigenvalue**3-eigenvalue**2*(d11+d22+d33)+eigenvalue*(d11*d22-d12*d21+d11*d21-d13*d31+d22*d33-d23*d32)+(-d13*d21*d32-d12*d23*d31-d11*d22*d33+d13*d22*d31+d12*d21*d33+d11*d23*d32)
    if eqn >0 and eqn<1:
        sol_eigen.append(eigenvalue)
i=0
for eig in sol_eigen:
    print('eigen'+i+'='+eig)
    i+=1

