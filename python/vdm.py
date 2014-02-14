from __future__ import print_function
from pylab import *

'''Vehicle Dynamics Model by Rakha et al.'''

#define physical constants:
g=9.8      #Earth gravitational acceleration, m/s^2

#conversion factors
ms2kmh=3.6
kmh2ms=1/ms2kmh

#data for diesel fuel
diesel_mass_density=0.754 # kg/L
diesel_energy_density=35.86e6 # J/L
diesel_engine_efficiency=0.5 # dimensionless

#tractive effort section
#definitions

Pmax=6.83e5#engine power spec., W
#V=20      #speed, km/h
tiny=1e-3  #softening factor to avoid divide by zero
eta=0.9   #transmission efficiency, dimensionless, [0.89-0.94] - 0.1 for accessories
Mt=67000   #vehicle mass - empty (kg)
Ml=0
Ml=98000   #load mass (kg)  --- max at 98000
M=Mt+Ml    #total mass (kg)
fta=0.64   #fraction on tractive axle, [0-1] dimensionless, 0.64 for dump truck
Mta= M*fta #mass on tractive axle, kg
mu=0.60    #coeffiecient of friction, dimensionless
           #value from hpwizard.com: offroad tyre on dry gravel road

#resistance forces section, N
#definitions
c1= 0.047      # factor to due air density
    # 0.5*density / 3.6^2 (to convert from km/h to m/s)
    # comes from Ra=(1/2)*(density of air)*Cd*A*V^2
    # density of air is a function of temperature, see 
    # http://en.wikipedia.org/wiki/Density_of_air
Cd= 0.9        # drag coefficient, dimensionless
H = 100        # altitude, m
Ch= 1-8.5e-5*H # altitdue coefficient, dimensionless
Area = 31.5    # area, m^2

# for the fancy version of rolling resistance
#Cr = 2.50 #rolling coefficient
#c2 = 0.1  #rolling resistance coefficient
#c3 = 10   #rolling resistance coefficient
#Rr= g*Cr*(c2*v+c3)*M/1e3 #rolling resistance (friction)

friction=0.0   # the way they do it in mining
grade=0.05      # slope


def minmax(x,a,b): #utility
    return max(a,min(x,b))

def power(x): # set engine power, x in [0,1]
    return Pmax*minmax(x,0,1)

def acceleration(v,x):
    '''calculates acceleration for a truck moving at speed V'''
    #input: v, km/h, velocity
    #input: x, [0,1], throttle

    #tractive effort, N
    Ft=ms2kmh*eta*power(x)/(v+tiny)
    Fmax=g*Mta*mu
    F=min(Ft,Fmax)

    #resistance forces
    Ra= c1*Cd*Ch*Area*v*v    #air drag
    Rg = g*M*(grade+friction)

    #total resistance
    R= Ra+Rg

    return (F-R)/M #acceleration, m/s^2


def drive_truck(tmax=3600,v_target=-1):
    #input 
    # tmax = time to stop (s)
    # v_target = target velocity, if > 0 then adjust throttle

    use_throttle = v_target > 0

    print(" v_target = {0:4.1f}".format(v_target))

    t =0.   # time (s)
    dt=.1   # time step (s)
    v=0.    # velocity (km/h)
    dist=0. # distance (km)
    f=1135. # fuel (L)
    x=1.0   # throttle [0,1] dimensionless
    #initialize storage arrays (which are actually lists)
    fuel=[f]
    velocity=[v]
    distance=[dist] #km
    time=[t]
    throttle=[x]
    q=0.2 #sets width of interval over which throttle is adjusted

    #simple Euler's method ODE integration
    while t <= tmax:
        a=acceleration(v,x) # input v in km/h, output a in m/s^2
        dist+=1e-3*dt*(kmh2ms*v + 0.5*a*dt) # update distance (km)
        v+=a*dt*ms2kmh      # update velocity, convert from m/s to km/h
        t+=dt

        #fuel burned as energy balance 
        output=power(x)*dt #W*s=J
        input_per_L=diesel_engine_efficiency*diesel_energy_density
        fuel_burned = output/input_per_L # W / (W/L) = L
        f-=fuel_burned
        
        velocity.append(v)
        distance.append(dist)
        time.append(t)
        fuel.append(f)
        throttle.append(x)

        # adjust throttle by comparing v with v_target
        # throttle control model is Fermi-Dirac fxn with
        # parameter q to control the width of the interval
        if use_throttle: x=1./(exp((v-v_target)/q)+1.)
        #if use_throttle:
        #    if v > v_target:
        #        x*=0.8
        #    elif v < v_target:
        #        x*=1.2
        #    x=minmax(x,0,1)

    print(" v_final = {0:4.1f}".format(velocity[-1]))
    return array(time),array(velocity),array(distance),array(fuel),array(throttle)

close("all")
fs=20



tmax=200 #(s)

grade=0.05
t1,v1,d1,f1,thr1=drive_truck(tmax)


grade=0.10
t2,v2,d2,f2,thr2=drive_truck(tmax)


grade=0.15
t3,v3,d3,f3,thr3=drive_truck(tmax)




#do one with a target velocity
#v_target=10. #kph
#t2,v2,d2,f2,thr2=drive_truck(tmax,v_target)

if False:
    xmin=0
    xmax=tmax

    #plot quantities vs. time
    fig1=figure(1,figsize=(10,10))
    fig1.subplots_adjust(bottom=0.08,top=0.96,right=0.96)
    subplot(311)
    plot(t1,v1,color="Blue",label="grade 5%")
    plot(t2,v2,color="Red",label="grade 10%")
    plot(t3,v3,color="Green",label="grade 15%")
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    ylabel("velocity (km/h)",fontsize=fs)
    legend(loc='lower right',fancybox=True,fontsize=fs)
    xlim(xmin,xmax)

    subplot(312)
    plot(t1,d1,color="Blue")
    plot(t2,d2,color="Red")
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    ylabel("distance (km)",fontsize=fs)
    xlim(xmin,xmax)

    #subplot(413)
    #plot(t1,f1,color="Blue")
    #plot(t2,f2,color="Red")
    #xticks(fontsize=fs)
    #yticks(fontsize=fs)
    #ylabel("fuel (L)",fontsize=fs)
    #xlim(xmin,xmax)

    subplot(313)
    plot(t1,thr1,color="Blue")
    plot(t2,thr2,color="Red")
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    xlabel("time (s)",fontsize=fs)
    ylabel("throttle",fontsize=fs)
    xlim(xmin,xmax)
    show()
