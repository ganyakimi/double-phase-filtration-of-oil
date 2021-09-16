import numpy as np
import matplotlib.pyplot as plt


def f(s):
    return (s*s/(s*s+0.1*(1-s)*(1-s)))

def f1(s):
    return (s*s) 

def f2(s):
    return (1-s)*(1-s)
def k(x):
    return 4

N = 150

M = 700

s=np.zeros((N+1,M+1),'float')
p=np.zeros((N+1,M+1),'float')
a = np.zeros((N+1),'float')
al = np.zeros((N+1),'float')
be = np.zeros((N+1),'float')
x=np.zeros((M+1),'float')

h = 0.01
tau = 0.002

#q, x

 

for i in range (N+1):
    s[i][0]=0.0
    x[i]=i*h

for j in range (M):
    s[0][j+1]=1.0
    al[0]= 1.0
    be[0]= 1.0

    for i in range (N):
        a[i+1]=k(x[i])*(f1((s[i][j]+s[i+1][j])/2)+1/10.*f2((s[i][j]+s[i+1][j])/2))
    for i in range (1,N):
        al[i]=a[i]*al[i-1]/(a[i+1]+al[i-1]*a[i])
        be[i]=a[i]*be[i-1]/(a[i+1]+al[i-1]*a[i])

    p[100][j+1]=0

    for i in reversed(range(0, N-1)):
        p[i][j+1]=(1-al[i])*p[i+1][j+1]+be[i]
    q=-a[99]*(p[99][j+1]-p[98][j+1])/h
    for i in range (1,N+1):
        s[i][j+1]=s[i][j]-tau/h*q*(f(s[i][j])-f(s[i-1][j]))

s1=np.zeros((N+1),'float')
s2=np.zeros((N+1),'float')
s3=np.zeros((N+1),'float')

p1=np.zeros((N+1),'float')
p2=np.zeros((N+1),'float')
p3=np.zeros((N+1),'float')
 
for i in range (N+1):
    print(s[i][j])
for i in range (N+1):
    print(p[i][j])

for j in range (N+1):
        s1[j]=s[j][200]
        s2[j]=s[j][400]
        s3[j]=s[j][650]
        p1[j]=p[j][200]
        p2[j]=p[j][400]
        p3[j]=p[j][650]

plt.figure(1)
plt.plot(s1)
plt.plot(s2)
plt.plot(s3)
plt.xlim(-5,100)
plt.xlabel('x')
plt.ylabel('s')
plt.legend(['t=200', 't=400', 't=650'])
plt.savefig('s_odnorod n+v.png')
plt.show()
plt.figure(2)
plt.plot(p1)
plt.plot(p2)
plt.plot(p3)
plt.legend(['t=200', 't=400', 't=650'])
plt.xlabel('x')
plt.ylabel('p')
plt.savefig('p_odnorod n+v.png')
plt.show()
plt.close('all')
