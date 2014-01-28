from pylab import *

'''Vehicle Dynamics Model by Rakha et al.'''

#define physical constants:
g=9.8      #Earth gravitational acceleration, m/s^2

#data for diesel fuel
diesel_mass_density=0.754 # kg/L
diesel_energy_density=35.86e6 # J/L
diesel_engine_efficiency=0.5 # dimensionless

#class Truck(object):
#''' a truck '''
#    def __init__(self):
#        self.power=power  #engine power (W)
#        self.trans_eff=tran_eff #transmission efficiency, dimensionless, [0-1]
#        self.eng_eff=eng_eff
#        self.vehicle_mass=vehicle_mass #unloaded mass of the vehicle (kg)
#        self.load_mass=load_mass #mass of load (kg)
#        self.fta=fta #fraction on tractive axle, dimensionless, [0-1]
#        self.drag=drag #drag coefficient, dimensionless
#        self.area=area #vehicle front surface area, m^2
#        self.position=position # [x,y,z]
#        self.fuel_tank_max=1135 # L (optional 1324 L)
#        self.fuel_tank=0


#tractive effort section
#definitions

P=683      #engine power, kW
#V=20       #speed, km/h
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


def acceleration(V):
    '''calculates acceleration for a truck moving at speed V'''
    #input: V, km/h

    #tractive effort, N
    Ft=3600*eta*P/(V+tiny)
    Fmax=g*Mta*mu
    F=min(Ft,Fmax)

    #resistance forces
    Ra= c1*Cd*Ch*Area*V**2      #air drag
    Rr= g*Cr*(c2*V+c3)*M/1e3 #rolling resistance (friction)
    Rg= g*M*grade            #grade resistance

    #total resistance force
    R = Ra + Rr + Rg

    return (F-R)/M #acceleration, m/s^2



close("all")
fs=20
if True:
    #test it:
    t =0. # time (s)
    dt=1. # time step (s)
    v=0.  # velocity (km/h)
    f=1135 # fuel (L)
    fuel=[f]
    velocity=[v]
    time=[t]
    #simple Euler's method ODE integration for one hour.
    while t < 3600:
        a=acceleration(v)   # input v in km/h, output a in m/s^2
        #if a < 1e-4: break # stop when a gets too small.
        v+=a*dt/3.6         # update velocity, convert from m/s to km/h
        t+=dt
        
        #fuel burned as energy balance 
        output=P*1e3*dt #W*s=J
        input_per_L=diesel_engine_efficiency*diesel_energy_density
        fuel_burned = output/input_per_L # W / (W/L) = L
        f-=fuel_burned
        
        velocity.append(v)
        time.append(t)
        fuel.append(f)
        
    time=array(time)
    velocity=array(velocity)
    fuel=array(fuel)

    #plot velocity vs. time
    figure(1,figsize=(8,8))
    plot(time,velocity)
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    xlabel("time (s)",fontsize=fs)
    ylabel("velocity (km/h)",fontsize=fs)
    
    figure(2,figsize=(8,8))
    plot(time,fuel/fuel[0])
    xticks(fontsize=fs)
    yticks(fontsize=fs)
    xlabel("time (s)",fontsize=fs)
    ylabel("fuel (L)",fontsize=fs)
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
