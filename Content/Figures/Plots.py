import matplotlib.pyplot as plt

T = []
RHO0 = []
P0 = []
P1 = []

indata=open('out.txt','r')
for aline in indata:
	t,rho0,p0,p1,=aline.strip().split()
	T.append(t)
	RHO0.append(rho0)
	P0.append(p0)
	P1.append(p1)
indata.close()

plt.plot(T,P0, label = 'p0')
plt.plot(T,P1, label = 'p1')
plt.xlabel('Time', fontsize=14, fontweight='bold')
plt.ylabel('Pressure', fontsize=14, fontweight='bold')
plt.legend(loc='upper left')
plt.legend()
plt.show()