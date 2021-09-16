import numpy as np
import matplotlib.pyplot as plt


def mu(c):
	return 20 /(1 + 5*c)

def f1(s):
	return (s*s)
 
def f2(s):
	return (1 - s) *(1-s)

def f(s,c):
	return s*s/ (s*s + mu(c)*(1-s)*(1-s))

def k(x):
	if(x<0.5):
		return 1
	else:
		return 5

N = 100
M = 900

s=np.zeros((N+1,M+1),'float')
p=np.zeros((N+1,M+1),'float')
a = np.zeros((N+1),'float')
al = np.zeros((N+1),'float')
be = np.zeros((N+1),'float')
x=np.zeros((M+1),'float')
c=np.zeros((N+1,M+1),'float')

h = 0.01
tau = 0.001
m=0.2

#q, x

for i in range (N+1):
	c[i][0] = 0.0
	s[i][0]=0.0
	x[i]=i*h

for j in range (M):
	s[0][j+1]=1.0
	al[0]= 1.0
	be[0]= 1.0
	c[0][j]=1.

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
		s[i][j+1]=s[i][j]-tau/h*q*(f(s[i][j],c[i][j])-f(s[i-1][j],c[i-1][j]))
		ct = s[i][j]*c[i][j] + 0.1*c[i][j] -q*(f(s[i][j], c[i][j])*c[i][j] - f(s[i-1][j], c[i-1][j])*c[i-1][j])*tau/h
		c[i][j+1] = ct/(s[i][j+1] + 0.1)
s1=np.zeros((N+1),'float')
s2=np.zeros((N+1),'float')
s3=np.zeros((N+1),'float')

p1=np.zeros((N+1),'float')
p2=np.zeros((N+1),'float')
p3=np.zeros((N+1),'float')

c1=np.zeros((N+1),'float')
c2=np.zeros((N+1),'float')
c3=np.zeros((N+1),'float')

for i in range (N+1):
        print(s[i][j])

for i in range (N+1):
	print(p[i][j])

for i in range (N+1):
	print(c[i][j])

for j in range (N+1):
        s1[j]=s[j][300]
        s2[j]=s[j][600]
        s3[j]=s[j][900]
        p1[j]=p[j][300]
        p2[j]=p[j][600]
        p3[j]=p[j][900]
        c1[j]=c[j][300]
        c2[j]=c[j][600]
        c3[j]=c[j][899]

plt.figure(1)
plt.plot(s1)
plt.plot(s2)
plt.plot(s3)
plt.ylim(0.8,1)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('s')
plt.savefig('s2.png')
plt.show()
plt.figure(2)
plt.plot(p1)
plt.plot(p2)
plt.plot(p3)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('p')
plt.savefig('p2.png')
plt.show()
plt.figure(3)
plt.plot(c1)
plt.plot(c2)
plt.plot(c3)
plt.legend(['t=300','t=600','t=900'])
plt.xlabel('x')
plt.ylabel('c')
plt.savefig('c2.png')
plt.show()
plt.close('all')
