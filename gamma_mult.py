import numpy as np
import matplotlib.pyplot as plt
import gmpy2 as g2
def diagnorm(h):
  return np.sum(np.abs(h[:,1])**2)
def totalnorm(h):
  return np.sum(np.abs(h)**2)
def bincount(k):
  #return bin(k).count("1")
  return g2.popcount(k)
def mult_fac(p,q,r,s):
  ans = bincount((s^q)&(p&r))
  ans += bincount((p^r)&(q&s))
  ans *= 2
  ans += bincount(r&q)
  ans -= bincount(p&s)
  return 1.j**ans
def unop(theta,s1,s2,hnew,h,n):
  hnew = h*(np.cos(theta)**2)
  for x in range(n):
    for y in range(n):
      hnew[x,y]+=h[x,y]*(np.sin(theta)**2)*mult_fac(s1,s2,s1^x,s2^y) \
                *mult_fac(x,y,s1,s2)
      hnew[x^s1,y^s2] += 0.5j*h[x,y]*np.sin(2*theta) \
                        *(mult_fac(s1,s2,x,y)-mult_fac(x,y,s1,s2))
  return hnew
n = 2**3
h = np.array(np.random.rand(n,n), dtype=complex)
h0=h+0.0
#h = np.array([[0,1,1,0],[0,0,0,0],[0,0,1,0],[0,1,0,0]],dtype=complex)
hnew = np.zeros((n,n),dtype=complex)
nth = 2**6
theta=np.linspace(0,2*np.pi,nth)
norms = np.zeros(nth,dtype=np.double)
def findmaxtheta(theta,s1,s2,h,hnew,nth,n):
  for i in range(nth):
    hnew=unop(theta[i],s1,s2,hnew,h,n)
    norms[i]=diagnorm(hnew)
  ans=0
  if np.max(abs(np.diff(norms)))>0.01 :
    ans=norms.argmax()
  else:
    ans = np.random.randint(0,nth-1)
  return ans
m = 11
pdata=[]
for j in range(m):
  for x in range(n):
    for y in range(n):
      a = findmaxtheta(theta,x,y,h,hnew,nth,n)
      h=unop(theta[a],x,y,hnew,h,n)
      pdata.append(diagnorm(h)/totalnorm(h))
plt.plot(pdata)
plt.ylabel('some numbers')
plt.show()


























