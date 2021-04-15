import math


print('insert d11[m]')
d11=float(input())
print('insert d12[m]')
d12=float(input())
print('insert d21[m]')
d21=float(input())
print('insert d22[m]')
d22=float(input())

print('w^4+'+str((-d11-d22))+'*w^2+'+str(d11*d22-d12*d21))

b=d11+d22
a=1
c=d11*d22-d12*d21


d = (b ** 2) - (4 * a * c)

# find two solutions
sol1 = (-b - math.sqrt(d)) / (2 * a)
sol2 = (-b + math.sqrt(d)) / (2 * a)

sol1 = math.sqrt(math.fabs(sol1))
sol2 = math.sqrt(math.fabs(sol2))
print('The solution are {0} rad/s and {1} rad/s'.format(sol1, sol2))


a_out="{:e}".format(a)
b_out="{:e}".format(b)
c_out="{:e}".format(c)

print("a: "+str(a_out))
print("b: "+str(b_out))
print("c: "+str(c_out))