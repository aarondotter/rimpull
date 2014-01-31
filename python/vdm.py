from __future__ import print_function
from pylab import *

'''Vehicle Dynamics Model by Rakha et al.'''

#define physical constants:
g=9.8      #Earth gravitational acceleration, m/s^2

#data for diesel fuel
diesel_mass_density=0.754 # kg/L
diesel_energy_density=35.86e6 # J/L
diesel_engine_efficiency=0.5 # dimensionless

#tractive effort section
#definitions

Pmax=6.83e5#engine power spec., W
#V=20      #speed, km/h
tiny=1e-3  #softening factor to avoid divide by zero
eta=0.8    #transmission efficiency, dimensionless, [0.89-0.94] - 0.1 for accessories
M=165000   #vehicle mass (kg)
fta=0.64   #fraction on tractive axle, [0-1] dimensionless, 0.64 for dump truck
Mta= M*fta #mass on tractive axle, kg
mu=0.3      #coeffiecient of friction, dimensionless


#resistance forces section, N
#definitions
c1= 0.047      # factor to due air density
    # 1/2*density / 3.6^2 (to convert from km/h to m/s)
    # comes from Ra=(1/2)*(density of air)*Cd*A*V^2
    # density of air is a function of temperature, see 
    # http://en.wikipedia.org/wiki/Density_of_air
Cd= 0.7        # drag coefficient, dimensionless
H = 100        # altitude, m
Ch= 1-8.5e-5*H # altitdue coefficient, dimensionless
Area = 31.5    # area, m^2

Cr = 2.50 #rolling coefficient
c2 = 0.1  #rolling resistance coefficient
c3 = 10   #rolling resistance coefficient

grade=0. # vertical / horizontal

def minmax(x,a,b): #utility
    return max(a,min(x,b))

def power(x): # set engine power, x in [0,1]
    return Pmax*minmax(x,0,1)

def acceleration(v,x):
    '''calculates acceleration for a truck moving at speed V'''
    #input: v, km/h, velocity
    #input: x, [0,1], throttle

    #tractive effort, N
    Ft=3.6*eta*power(x)/(v+tiny)
    Fmax=g*Mta*mu
    F=min(Ft,Fmax)

    #resistance forces
    Ra= c1*Cd*Ch*Area*v*v    #air drag
    Rr= g*Cr*(c2*v+c3)*M/1e3 #rolling resistance (friction)
    Rg= g*M*grade            #grade resistance

    #total resistance force
    R = Ra + Rr + Rg

    return (F-R)/M #acceleration, m/s^2


def drive_truck(tmax=3600,v_target=-1):
    #input 
    # tmax = time to stop (s)
    # v_target = target velocity, if > 0 then adjust throttle

    use_throttle = v_target > 0

    t =0.  # time (s)
    dt=1.  # time step (s)
    v=0.   # velocity (km/h)
    f=1135 # fuel (L)
    x=1.0  # throttle [0,1] dimensionless
    v_target = 30. # target velocity (km/h)
    fuel=[f]
    velocity=[v]
    time=[t]
    throttle=[x]
    q=1e-1 #sets width of interval over which throttle is adjusted

    #simple Euler's method ODE integration for one hour.
    while t <= tmax:
        a=acceleration(v,x) # input v in km/h, output a in m/s^2
        #if a < 1e-4: break # stop when a gets too small.
        v+=a*dt/3.6         # update velocity, convert from m/s to km/h
        t+=dt

        #fuel burned as energy balance 
        output=power(x)*dt #W*s=J
        input_per_L=diesel_engine_efficiency*diesel_energy_density
        fuel_burned = output/input_per_L # W / (W/L) = L
        f-=fuel_burned
        
        velocity.append(v)
        time.append(t)
        fuel.append(f)
        throttle.append(x)

        # adjust throttle by comparing v with v_target
        # throttle control model is Fermi-Dirac with
        # parameter q to control the width of the interval
        if use_throttle: x=1./(exp((v-v_target)/q)+1.)

    return array(time),array(velocity),array(fuel),array(throttle)

close("all")
fs=20

#do one full throttle
tmax=1200 #(s)
t1,v1,f1,thr1=drive_truck(tmax)

#do one with a target velocity
v_target=30 #kph
t2,v2,f2,thr2=drive_truck(tmax,v_target)

xmin=0
xmax=tmax
#plot velocity vs. time
fig1=figure(1,figsize=(10,10))
fig1.subplots_adjust(top=0.96,right=0.96)
subplot(311)
plot(t1,v1,color="Blue",label="full")
plot(t2,v2,color="Red",label="adjust")
xticks(fontsize=fs)
yticks(fontsize=fs)
#xlabel("time (s)",fontsize=fs)
ylabel("velocity (km/h)",fontsize=fs)
legend(loc='lower right',fontsize=fs)
xlim(xmin,xmax)

#figure(2,figsize=(8,8))
subplot(312)
plot(t1,f1,color="Blue")
plot(t2,f2,color="Red")
xticks(fontsize=fs)
yticks(fontsize=fs)
#xlabel("time (s)",fontsize=fs)
ylabel("fuel (L)",fontsize=fs)
xlim(xmin,xmax)

#figure(3,figsize=(8,8))
subplot(313)
plot(t1,thr1,color="Blue")
plot(t2,thr2,color="Red")
xticks(fontsize=fs)
yticks(fontsize=fs)
xlabel("time (s)",fontsize=fs)
ylabel("throttle",fontsize=fs)
xlim(xmin,xmax)
ylim(0.8,1.01)
show()

if False:
    masses=linspace(75000,165000,10)
    grades=linspace(0,0.30,7)
    
    results=[]
    for M in masses:
        res=[]
        for grade in grades:
            t=0
            dt=1
            v=0
            while t < 3600:
                a=acceleration(v)
                if a < 1e-5: break
                v+=a*dt/3.6
                t+=dt
            #print M, grade, v
            res.append(v)
        results.append(res)
    results=array(results)
    
    figure(2,figsize=(8,8))
    for i,grade in enumerate(grades):
        plot(masses*1e-3,results[:,i],label="{0:4.2f}".format(grade))
    xlim(70,170)
    ylim(0,70)
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    xlabel("gross vehicle weight x 1000 (kg)",fontsize=fs)
    ylabel("maximum velocity (km/h)",fontsize=fs)
    legend(fancybox=True,title="Grade",fontsize=fs-2,
           borderpad=0.2,handletextpad=0.2,labelspacing=0.2)
    show()
