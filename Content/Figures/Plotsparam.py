import matplotlib.pyplot as plt

P1 = []
P2 = []
LNR = []
count = 1
indata=open('output.txt','r')
for aline in indata:
	LNR.append(count) 
	param1, param2 =aline.strip().split()
	P1.append(param1)
	P2.append(param2)
	count = count + 1
indata.close()


fig, ax1 = plt.subplots()
ax1.plot(LNR,P1, label = 'P1', color = 'b')
ax1.set_ylabel('Parameter 1', color = 'b')
ax1.set_xlabel('Iternation Number', color = 'g')
ax2 = ax1.twinx()
ax2.plot(LNR,P2, label = 'P2', color = 'r')
ax2.set_ylabel('Parameter 2', color = 'r')
plt.show()