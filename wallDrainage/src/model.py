#!/usr/bin/env python
# TODO: fix diffusion equation - dimensionless position
from __future__ import division
from fipy import Grid1D,CellVariable,Variable,TransientTerm,DiffusionTerm,\
    Viewer
from scipy.integrate import ode
from numpy import zeros,exp,sqrt
import datetime
neval=[]
ts = datetime.datetime.now()
# create mesh
nx = 200
dx = 1.e-6
mesh = Grid1D(nx=nx, dx=dx)
# define variable
hf = CellVariable(name="concentration",mesh=mesh,hasOld=1)
# set initial condition
rd=100e-6
hi=5e-6
s=1/sqrt(3.0)
dstr=0.2
rs=rd*sqrt((1+s**2)/s**2)*dstr
for i in range(nx):
    hf[i]=(rs+hi)-sqrt(rs**2-(r(i)-(1-dstr)*rd)**2)
conc.setValue(conc0)
XOH0=0.
XW0=0.
temp0=304
R0=1e-6
p0=zeros(2)
p0[0]=0.
p0[1]=0.
# diffusion coefficient
def diff_model(D0,conc):
    return D0 * (1 + conc)
D0 = 1.e-12
# equation
eq = TransientTerm(var=conc) == \
    DiffusionTerm(coeff=diff_model(D0,conc),var=conc)
# show result when we run script directly
# if __name__ == '__main__':
#     viewer_conc = Viewer(vars=conc,datamin=0.,datamax=2.)
    # viewer_temp = Viewer(vars=temp,datamin=200.,datamax=400.)
# timestep
def odesystem(t,y):
    xOHeq=0
    xWeq=1
    teq=2
    req=3
    fpeq=4
    lpeq=5
    AOH=1965
    EOH=5.487e4
    dHOH=-7.49e4
    AW=1.385e3
    EW=3.266e4
    dHW=-8.6e4
    NCO0=3765
    OH0=3765
    W0=0
    rhop=1100.
    cp=1800.
    Rg=8.314
    pamb=1.01e5
    sigma=25e-3
    D=zeros(2)
    D[0]=2.4e-12
    D[1]=4.4e-12
    Mbl=zeros(2)
    Mbl[0]=137.37e-3
    Mbl[1]=44e-3
    KH=zeros(2)
    KH[0]=rhop/Mbl[0]/pamb*(1e-7+4.2934*\
        exp(-(y[teq]-203.3556)**2/(2*40.016**2)))
    KH[1]=1.1e-4
    Aeta=4.1e-8
    Eeta=38.3e3
    Cg=0.85
    AA=4.e0
    B=-2.e0
    eta=Aeta*exp(Eeta/(Rg*y[teq]))*(Cg/(Cg-y[xOHeq]))**(AA+B*y[xOHeq])
    pair0=(pamb+2*sigma/y[req])
    pair=pair0*R0**3/y[req]**3*y[teq]/temp0
    ydot=zeros(6)
    ydot[xOHeq]=AOH*exp(-EOH/Rg/y[teq])*(1-y[xOHeq])*\
        (NCO0-2*W0*y[xWeq]-OH0*y[xOHeq])
    ydot[xWeq] = AW*exp(-EW/Rg/y[teq])*(1-y[xWeq])
    ydot[teq]=-dHOH*OH0/(rhop*cp)*ydot[xOHeq]-dHW*W0/(rhop*cp)*ydot[xWeq]
    ydot[req] = (sum(y[fpeq:lpeq]) + pair - pamb - \
        2*sigma/y[req])*y[req]/(4*eta)
    for i,j in enumerate(range(fpeq,lpeq+1)):
        ydot[j]=-3*y[j]*ydot[req]/y[req] + y[j]/y[teq]*ydot[teq] + \
            9*Rg*y[teq]*D[i]*y[req]*(conc[0]-KH[i]*y[j])/(dx/2)
    neval.append(1)
    return ydot

# define ODE integrator
solver=ode(odesystem).set_integrator('vode',method='bdf')
solver.set_initial_value([XOH0,XW0,temp0,R0,p0[0],p0[1]], 0.)
# boundary conditions
temp = Variable()
conc.constrain(1100/137.37e-3/1.01e5*(1e-7+4.2934*\
    exp(-(temp-203.3556)**2/(2*40.016**2)))*solver.y[4], mesh.facesLeft)
conc.faceGrad.constrain([0.], mesh.facesRight)
dt = 2.e0
steps = 50
# define number of sweeps
sweeps=20
res=1e8
print(solver.t,solver.y)
for step in range(steps):
    while 1 in neval:
        neval.remove(1)
    solver.integrate(solver.t+dt)
    # print len(neval)
    print(solver.t,solver.y)
    conc.updateOld()
    temp.setValue(solver.y[2])
    for sweep in range(sweeps):
        res = eq.sweep(dt=dt)
        if res < 1e-6:
            # print 'exiting after {0} sweeps'.format(sweep)
            break
    # if __name__ == '__main__':
    #     viewer_conc.plot()
        # viewer_temp.plot()

tf = datetime.datetime.now()
te = tf - ts
print "Runtime: ",te
if __name__ == '__main__':
    raw_input("Residual = %f. Press <return> to proceed..." % (abs(res)))
