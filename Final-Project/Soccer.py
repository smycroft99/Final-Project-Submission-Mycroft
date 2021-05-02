#python Documents\PHY2200\Final-Project\Soccer.py

#This program is intended to accurately model the trajectory of one or more soccer ball kicked from
#the penalty spot on a soccer field at a goal. The program plots each modeled ball's trajectory as
#it approaches the goal and colors the line of each of these balls dependent on it's initial spin
#conditions and the user's requested modeling parameters. For example, the user may ask to model
#100 soccer balls shot at the upper left corner of the goal and color them depending on the magnitude of their
#top/backspin.

#Import functions to use during program
import numpy as np
import math
import matplotlib.pyplot as plt
import random as rand
import ode
import ipywidgets as widgets

from IPython.display import display
from vpython import *
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

#Function definitions
def mag(v):           #Calculate magnitude of an array
    return np.sqrt(np.dot(v,v))

def hat(v):           #Calculate unit vector of an array
    return v/mag(v)

def cross(v1,v2):     #Calculate cross product of two arrays
    return np.cross(v1,v2)

def getCd(v):         #Return the value of Cd depending on the speed of the ball
    a = 0.36
    b = 0.14
    c = 0.27
    vc = 34
    chi = (v - vc)/4
    if chi < 0:
        Cd = a + b/(1+np.exp(chi)) - c*np.exp(-chi**2)
    else:
        Cd = a + b/(1+np.exp(chi)) - c*np.exp(-chi**2/4)
    
    return Cd

def model_magnus_2(d, tn): #Return array of derivatives
                           #rate[3]-rate[5] give accelerations due to magnus, drag, and gravity
        
    #data
    x = d[0]
    y = d[1]
    z = d[2]
    vx = d[3]
    vy = d[4]
    vz = d[5]
    
    #derivatives
    rate = np.zeros(6) #derivatives
    rate[0] = vx
    rate[1] = vy
    rate[2] = vz
    
    #speed
    v = np.array([vx,vy,vz])
    vmag = mag(v)
    
    #calculate magnus force
    Cd=getCd(vmag)
    b2 = 1/2*Cd*rho*A
    
    S = r*omegamag/vmag
    CL = 0.62*S**0.7
    alpha = 1/2*CL*rho*A*r/S
    
    FM = alpha*cross(omega[i],v) #magnus force
    
    speed = np.sqrt(rate[0]**2+rate[1]**2+rate[2]**2)
    
    #calculate dv/dt due to all forces
    #          [Linear Drag]  [Magnus]  [velocity dependent drag]                         [grav]
    rate[3] = -b2*vmag*vx/m + FM[0]/m - (6*np.pi*mu*r*speed*vx+0.5*Cd*rho*A*speed*vx)/m
    rate[4] = -b2*vmag*vy/m + FM[1]/m - (6*np.pi*mu*r*speed*vy+0.5*Cd*rho*A*speed*vy)/m - g
    rate[5] = -b2*vmag*vz/m + FM[2]/m - (6*np.pi*mu*r*speed*vz+0.5*Cd*rho*A*speed*vz)/m
    
    #return derivative values
    return rate

def run_magnus_2(data): # run simulation
    global b2, alpha, i, bounceDamp # need to change the value of b2 and alpha
    
    v = data[3:6]
    vmag = mag(v)

    #time
    t = 0
    h= 0.01
    Nsteps = int(10/h)

    #store trajectory for plotting or animation
    traj = np.zeros((Nsteps, 4)) #store t, x, y, z for plotting
    traj[0,:] = np.array([t, x0, y0, z0]) #store initial time and position

    for n in range(0,Nsteps-1):

        #update data
        data = ode.RK4(model_magnus_2, data, t, h )
        
        #Bounce off ground and lose energy
        if(data[1]<0):
            data[1] = -data[1]
            data[3] = data[3]*bounceDamp
            data[4] = -data[4]*bounceDamp
            data[5] = data[5]*bounceDamp

        #update t
        t = t + h

        #store data for plotting
        traj[n+1,:] = np.array([t, data[0], data[1], data[2]])
        if(data[0]>dist):
            break
    
    return traj

#Take in user input to determine:
#  a.)What color modeling will be used (darker colors reflect a higher magnitude of spin)
#Will not drop out until all answers are valid
while(1):
    print("What spin would you like to see?")
    print("[1] - Spin in direction of motion")
    print("[2] - Side spin")
    print("[3] - Top spin")
    print("[4] - All spin")
    UserInput = input()
    case = int(UserInput)-1
    if(case!=0 and case!=1 and case!=2 and case!=3):
        print("That is not a valid input, please enter 1, 2, 3, or 4.")
        print("")
    else:
        break
    
print("Magnitude of ball's spin in rotations per minute (600rpm is a good approximate value)")
UserInput = input()
spinrate = int(UserInput)
if(spinrate==0):
    spinrate = 1e-10

while(1):    
    print("Magnitude of ball's velocity immediately after kick in miles per hour (70mph is a good approximate value)")
    UserInput = input()
    v_pen = float(UserInput)
    if (v_pen<=0):
        print("Velocity must be greater than 0mph")
    else:
        break
while(1):
    print("How much loft would you like to give the ball? (Loft is the vertical angle relative to the ground where 0=flat along the ground and 90=directly vertical.) A standard shot at 70mph will hit the cross bar at an angle of ~17-18 degrees")
    print("WARNING: A higher loft angle will take longer to simulate due to the iteration setup")
    UserInput = input()
    theta = float(UserInput)
    if (theta<0 or theta>=90):
        print("Angle must be between 0 and 90 degrees to aim at the goal")
    else:
        theta = theta * np.pi/180
        break

while(1):
    print("How much range would you like to give the ball? (Range is the horizontal angle relative to the center of the goal where 0=center of the goal, -90=directly to your left, and 90=directly to your right.) With no spin, an angle of +/-18.4 degress will hit the goal post")
    UserInput = input()
    phi = float(UserInput)
    if (phi<-90 or phi>=90):
        print("Angle must be between -90 and 90 degrees to aim at the goal (0 for center of goal)")
    else:
        phi = phi * np.pi/180
        break
        
while(1):
    print("How many balls would you like to simulate? (Note: More shots will take longer to run. 100 shots takes ~25 sec to run)")
    UserInput = input()
    Nshots = float(UserInput)
    if(Nshots<=0 or Nshots%1!=0):
        print("Please enter a positive integer number of shots to simulate")
    else:
        Nshots=int(Nshots)
        break
                
while(1):
    print("Would you like to display lines that miss the goal?")
    print("[1] - Display all trajectories in full color")
    print("[2] - Display all trajectories that hit the goal in color, display trajectories that miss in black")
    print("[3] - Display all trajectories that hit the goal in color, display trajectories that miss at 2% opacity in color")
    UserInput = input()
    OnGoal = int(UserInput)
    if (OnGoal == 1 or OnGoal == 2 or OnGoal == 3):
        break
    else:
        print("Please enter Y for Yes or N for No")
       
#I wonder what happens if you enter the digits of pi one at a time as answers to the starting questions?       
if(case==2 and spinrate==1 and v_pen==4 and theta==1*(np.pi/180) and phi==5*(np.pi/180) and Nshots==9 and OnGoal==2):
    print("                           .-.")
    print("                          |_:_|")
    print("                         /(_Y_)\ ")
    print("    .                   ( \/M\/ )")
    print("     '.               _.'-/'-'\-'._")
    print("       ':           _/.--'[[[[]'--.\_")
    print("         ':        /_'  : |::\"| :  '.\ ")
    print("           ':     //   ./ |oUU| \.'  :\ ")
    print("             ':  _:'..' \_|___|_/ :   :|")
    print("               ':.  .'  |_[___]_|  :.':\ ")
    print("                [::\ |  :  | |  :   ; : \ ")
    print("                 '-'   \/'.| |.' \  .;.' |")
    print("                 |\_    \  '-'   :       |")
    print("                 |  \    \ .:    :   |   |")
    print("                 |   \    | '.   :    \  |")
    print("                 /       \   :. .;       |")
    print("                /     |   |  :__/     :  \\")
    print("               |  |   |    \:   | \   |   ||")
    print("              /    \  : :  |:   /  |__|   /|")
    print("              |     : : :_/_|  /'._\  '--|_\ ")
    print("              /___.-/_|-'   \  \ ")
    print("                  '-'")
    print("        COME TO THE DARK SIDE, WE HAVE EASTER EGGS")
    print("")
        

print("Running simulation, please wait...")

    
        #print("")
        #if(case==0):
        #    print("Generating a plot of %i soccer balls with randomized spin of magnitude %i rpm. The plot will depict shots from white to red depending on the magnitude of the absolute value of its parallel spin component." % (Nshots, spinrate))
        #if(case==1):
        #    print("Generating a plot of %i soccer balls with randomized spin of magnitude %i rpm. The plot will depict shots from white to green depending on the magnitude of the absolute value of its side spin component." % (Nshots, spinrate))
        #if(case==2):
        #    print("Generating a plot of %i soccer balls with randomized spin of magnitude %i rpm. The plot will depict shots from white to blue depending on the magnitude of the absolute value of its top/back spin component." % (Nshots, spinrate))
        #if(case==3):
        #    print("Generating a plot of %i soccer balls with randomized spin of magnitude %i rpm. The plot will depict shots with rgb values corresponding to the magnitude of the spin components where parallel=red, side spin=green, and top/back spin=blue." % (Nshots, spinrate))
        #print("")
    
L = 36 #ft                      #Distance from penalty spot to goal in ft
crossbar = 8 #ft                #Height of goal cross bar in ft

g = 9.8                         #N/kg
rho = 1.2                       #kg/m^3
mu = 1.8e-5                     #kg/m/s

r = 22e-2/2                     #220 mm diameter
A = np.pi*r**2                  #cross-sectional area
m = 0.4536                      #kg

dist = L/3.281                  #Convert distance to goal line to meters
CB = crossbar/3.281             #Convert crossbar height to meters
v_0 = v_pen/2.237               #Convert initial speed to m/s
omegamag = spinrate*np.pi/30    #Convert spin rate to radians per second
bounceDamp = 0.9                #Decimal value between 0 and 1 that represents energy conserved after
                                #ball bounces off ground

Nsteps = 1000                   #Number of steps each ball will take (1000 is good baseline)
t = 0                           #time
h = 0.01                        #dt

########################################################################################################
#Nshots = 5                      #How many balls will be simulated
#theta = math.atan(CB/dist)      #Define launch angle to hit center of goal (no gravity, drag, or magnus)
#theta = 0.3                     #Hardcode a value of theta
#phi = 0                         #Horizontal angle of shot (phi=0 for center of goal)
########################################################################################################

#Populate an array with [Nshots] random normalized angular velocity vectors of magnitude equal to [omegamag]
omega = []
for i in range(0,Nshots):    
    omX = rand.randint(-100,100)
    omY = rand.randint(-100,100)
    omZ = rand.randint(-100,100)
    omMag = np.sqrt(omX**2+omY**2+omZ**2)
    
    omX = omX/omMag*omegamag
    omY = omY/omMag*omegamag
    omZ = omZ/omMag*omegamag
    omMag = np.sqrt(omX**2+omY**2+omZ**2)
    
    omega.append([omX, omY, omZ])

#Initial conditions
x0 = 0
y0 = 0
z0 = 0
vx0 = v_0*np.cos(phi)*np.cos(theta)
vy0 = v_0*np.cos(phi)*np.sin(theta)
vz0 = v_0*np.sin(phi)

d = ([x0, y0, z0, vx0, vy0, vz0])
traj = []

for i in range(0,Nshots):
    temp = run_magnus_2(d)       #Simulate model for one ball with ith [omega] vector
     
    temp = temp[temp[:,2]>=0]    #get trajectory data with y>=0
    temp = temp[temp[:,1]<=dist] #get trajectory data with x<=distance to goalline
    traj.append(temp)

#Draw goal on 3D plot
GPx = np.array([[-3/2*CB, -3/2*CB], [3/2*CB, 3/2*CB]])
GPy = np.array([[dist, dist], [dist, dist]])
GPz = np.array([[0, CB], [0, CB]])
Fx  = np.array([[-CB*3, -CB*3], [CB*3, CB*3]])
Fy  = np.array([[-dist/4, dist], [-dist/4, dist]])
Fz  = np.array([[0, 0], [0, 0]])

#Produce array of hex-values corresponding to strength of desired omega axis
c2 = []
for j in range(Nshots):  
    Cout = '#'    
    ls=[k for k, e in enumerate(traj[j]) if e[2] != 0]
    last=traj[j][ls[-1]]
    print(last)
    for i in range(3):
        if(OnGoal==1 or OnGoal==3 or (OnGoal==2 and last[2]<CB and last[2]>=0 and np.abs(last[3])<CB*3/2)):
            if (case != i and case>=0 and case<=2 and OnGoal == 3 and (last[2]>=CB or last[2]<0 or np.abs(last[3])>=CB*3/2)):
                Cout = Cout + (format(255-abs(int(255*(omega[j][case])/(2*omegamag))), '02x'))
            elif (case != i and case>=0 and case<=2):
                Cout = Cout + (format(255-abs(int(255*(omega[j][case])/omegamag)), '02x'))
            elif (case<0 or case>2):
                Cout = Cout + (format(255-abs(int(255*(omega[j][i])/omegamag)), '02x'))
            else:
                Cout = Cout + (format(abs(255), '02x'))
        else:
            Cout = Cout + (format(abs(0), '02x'))
    c2.append(Cout)
#print(c2)

#Plot all [Nshots] shots on a 3D model with appropriate color values
ax = plt.axes(projection='3d', label='null')
plt.title("Soccer Ball Trajectory with Magnus Effect and Drag")
for i in range(Nshots):
    ls=[k for k, e in enumerate(traj[i]) if e[2] != 0]
    last=traj[i][ls[-1]]
    graphing=traj[i]
    if(OnGoal==3 and (last[2]>=CB or last[2]<0 or np.abs(last[3])>=CB*3/2)):
        ax.plot3D(graphing[:,3], graphing[:,1], graphing[:,2], color = c2[i], alpha=0.02)
    else:
        ax.plot3D(graphing[:,3], graphing[:,1], graphing[:,2], color = c2[i])
        
ax.plot_surface(Fx, Fy, Fz, color="#20F020", alpha=0.2)
ax.plot_surface(GPx, GPy, GPz, alpha=0.2)
ax.set_xlim3d(-dist/2,dist/2)
ax.set_ylim3d(0,dist)
ax.set_zlim3d(0,dist/2)

#Spin the graph slowly
angle = 0
while (plt.fignum_exists(1)):
    #angle = angle + 1
    #ax.view_init(30, angle)
    #plt.draw()
    plt.pause(.002)