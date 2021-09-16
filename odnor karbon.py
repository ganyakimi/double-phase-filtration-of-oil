import numpy as np
import matplotlib.pyplot as plt
import math


def mu1(c1):
    return ((5.5*c1*c1) - (14.7*c1) + 14.6)

def mu2(c2):
    return ((0.55 - (0.3*c2))*c2 + 1.0)

def f1(s):
    return s**2
 
def f2(s):
    return (1-s)**2

def f(s,c1,c2):
    mu=mu2(c2)/mu1(c1)
    return  (f1(s)/(f1(s) + (mu*f2(s))))

def k(x):   #проницаемость
    return 1


N = 600
M = 1000

s=np.zeros((N+1,M+1),'float')
p=np.zeros((N+1,M+1),'float')
a = np.zeros((N+1),'float')
al = np.zeros((N+1),'float')
be = np.zeros((N+1),'float')
x=np.zeros((M+1),'float')
c1=np.zeros((N+1,M+1),'float')
c2=np.zeros((N+1,M+1),'float')

h = 1./N
tau = 0.001
m=0.2
K=2.0 #коэффициент распределения 

#q, x

for i in range (N+1):
    x[i]=i*h
    c1[i][0] = 0.0
    c2[i][0]=0.0
    s[i][0]=0.0
    

for j in range (M):
    s[0][j+1]=1.0
    al[0]= 1.0
    be[0]= 1.0
    c1[0][j]=1.
    c2[0][j]=c1[0][j]/K

    for i in range (N):
        a[i+1]=k(x[i])*(f1((s[i][j]+s[i+1][j])/2)+1/10.*f2((s[i][j]+s[i+1][j])/2))

    for i in range (1,N):
        al[i]=a[i]*al[i-1]/(a[i+1]+al[i-1]*a[i])
        be[i]=a[i]*be[i-1]/(a[i+1]+al[i-1]*a[i])

    p[100][j+1]=0

    for i in reversed(range(0, N-1)):
        p[i][j+1]=(1-al[i])*p[i+1][j+1]+be[i]

    q=-a[99]*(p[99][j+1]-p[98][j+1])/h

    for i in range (1,N):
        c2[i][j] = c1[i][j] / K
        A = (K - ((K - 1) * f(s[i][j], c1[i][j], c2[i][j])))/(1.0 + (s[i][j]*(K - 1)))
        c1[i][j+1] =  c1[i][j] -(A *q* (c1[i][j] - c1[i-1][j]))
        s[i][j+1]=s[i][j]-tau*q/h*(f(s[i][j],c1[i][j],c2[i][j])-f(s[i-1][j],c1[i-1][j],c2[i-1][j]))

S1=np.zeros((N+1),'float')
S2=np.zeros((N+1),'float')
S3=np.zeros((N+1),'float')

P1=np.zeros((N+1),'float')
P2=np.zeros((N+1),'float')
P3=np.zeros((N+1),'float')

C11=np.zeros((N+1),'float')
C21=np.zeros((N+1),'float')
C31=np.zeros((N+1),'float')

C12=np.zeros((N+1),'float')
C22=np.zeros((N+1),'float')
C32=np.zeros((N+1),'float')

for i in range (N+1):
        print(s[i][j])

for i in range (N+1):
    print(p[i][j])

for i in range (N+1):
    print(c1[i][j])

for i in range (N+1):
    print(c2[i][j])

for j in range (N+1):
        S1[j]=s[j][300]
        S2[j]=s[j][600]
        S3[j]=s[j][900]
        P1[j]=p[j][300]
        P2[j]=p[j][600]
        P3[j]=p[j][900]
        C11[j]=c1[j][300]
        C21[j]=c1[j][600]
        C31[j]=c1[j][899]
        C12[j]=c2[j][300]
        C22[j]=c2[j][600]
        C32[j]=c2[j][899]

plt.figure(1)
plt.plot(S1)
plt.plot(S2)
plt.plot(S3)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('s')
plt.savefig('s3-1.png')
plt.show()
plt.figure(2)
plt.plot(P1)
plt.plot(P2)
plt.plot(P3)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('p')
plt.savefig('p3-1.png')
plt.show()
plt.figure(3)
plt.plot(C11)
plt.plot(C21)
plt.plot(C31)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('c1')
plt.savefig('c31-1.png')
plt.show()
plt.figure(4)
plt.plot(C12)
plt.plot(C22)
plt.plot(C32)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('c2')
plt.savefig('c32-1.png')
plt.show()
plt.close('all')

