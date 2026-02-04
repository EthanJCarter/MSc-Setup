#import modules
import numpy as np
import random

def convertSci(time,Ms,x,y,z,r,vx,vy,vz,v):

    #convert all values to scientific format for legibility
    timeconv = "{:.4e}".format(time)
    Msconv = "{:.4e}".format(Ms)
    xconv = "{:.4e}".format(x)
    yconv = "{:.4e}".format(y)
    zconv = "{:.4e}".format(z)
    rconv = "{:.4e}".format(r)
    vxconv = "{:.4e}".format(vx)
    vyconv = "{:.4e}".format(vy)
    vzconv = "{:.4e}".format(vz)
    vconv = "{:.4e}".format(v)

    return timeconv, Msconv, xconv, yconv, zconv, rconv, vxconv, vyconv, vzconv, vconv

def writeICs(inputfile, bodies, n_bodies):


    print("",len(bodies), "sinks", file = inputfile)

    print("id/pid/step/time(yr)/mass(Msun)/x/y/z/r(AU)/vx/vy/vz/v(km/s)", file = inputfile)

    for body in range (0,n_bodies):
        print(bodies[body], file = inputfile)

    #print("The input file is saved in this directory under the file name 'run_ic.dat'. ")



def makeIC_file(eccentricity, M_f, a_f, v_f, mp, p_a):
    
    inputfile = open("run_ic.dat", 'w')
    
    bodies, n_bodies = generate_bodies(eccentricity, M_f, a_f, v_f, mp, p_a)
    writeICs(inputfile,bodies,n_bodies)
    
    pass


def writeList(body_id, pid, step, time,Ms,x,y,z,r,vx,vy,vz,v, bodies):

    timeconv,Msconv,xconv,yconv,zconv,rconv,vxconv,vyconv,vzconv,vconv = convertSci(time,Ms,x,y,z,r,vx,vy,vz,v)

    list_entry = (" {} {} {} {} {} {} {} {} {} {} {} {} {}".format((body_id),pid,step,timeconv,Msconv,xconv,yconv,zconv,rconv,vxconv,vyconv,vzconv,vconv))

    #Appends temp list of current body properties to full list of bodies
    bodies.append(list_entry)

    #print(bodies)


def generate_bodies(eccentricity, M_f, a_f, v_f, mp, p_a):

    n_bodies = 3
    
    R4 = random.uniform(0,1)

    ta = 2*np.pi*R4 #True Anomaly, 2pi*R4, Radians

    for body_id in range(1, n_bodies+1):
        
        
        time = 0 #Yr
        pid = 0
        step = 0
        
        
        bodies = []
        
        if  body_id == 1:
            #Host Star
            Ms = 0.2 #Solar Masses
                    
            S_a = 0 #AU

            #Calculate initial positions based off of generated true anomaly
            S_x = 0 #Calculated from a and semi-major axis, AU
            S_y = 0 #As above, AU
            #S_z = 0
            S_r = 0

            #Calculate velocities based off of calculated positions and semi-major axis
            S_v = 0 #Total velocity (GM/R)^1/2, Km/s
            S_vx = 0 #-v*(y/a) #x-component, -vy/a, Km/s
            S_vy = 0 #v*(x/a) #y-component, -vx/a, Km/s
            #S_vz = 0

        elif body_id == 2:

            #Planet
            #Mp = random.uniform(0.00095,13*0.00095) #Solar Masses
            Mp = mp*0.00095
            #Mp = 1.2131e-02
                    
            #p_a =  random.randint(50,200) #AU

            a_metres = p_a*149597870700
            
                    
            p_r = (p_a*(1-eccentricity**2))/((1+(eccentricity*np.cos(ta))))
                                             

            r_metres = p_r*149597870700
            
            p_x = p_r*np.cos(ta) #Calculated from a and semi-major axis, AU
            p_y = p_r*np.sin(ta) #As above, AU
            #p_z = 0 
                    
            r_metres = 149597870700*p_r

            theta = np.pi/2 - ta
            
            #print(theta)

            p_v = np.sqrt((6.67408e-11*((Ms*2e30)+Mp))*((2/r_metres) - (1/a_metres)))/1000 #Total velocity (GM/R)^1/2, Km/s
            
            #print(p_v)
            
            p_vx = p_v * np.cos(theta) #/1000 #-p_v*(p_y/p_r) #x-component, -vy/a, Km/s
            p_vy = -p_v * np.sin(theta)  #/1000 #p_v*(p_x/p_r) #y-component, -vx/a, Km/s
            #p_vz = 0

    com_to_host,com_x,com_y,com_v,com_vx,com_vy = calc_com(ta, Ms, Mp, p_x, p_y, p_a, p_vx, p_vy, p_v, eccentricity)
    
    
    for body_id in range(1, n_bodies+1):
        
        if body_id == 1: #Host star in COM ref frame
            
            M = Ms
            
            a,x,y,z,r,vx,vy,vz,v = change_ref_frame(S_a,S_x,S_y,S_r,S_v,S_vx,S_vy,com_to_host,com_x,com_y,com_v,com_vx,com_vy)
          
        elif body_id == 2: #Planet in COM ref frame
            
            M = Mp
            
            a,x,y,z,r,vx,vy,vz,v = change_ref_frame(p_a,p_x,p_y,p_r,p_v,p_vx,p_vy,com_to_host,com_x,com_y,com_v,com_vx,com_vy)
            
            a_p = a
            v_p = v
        
        else:
        
        
            M = M_f
            a = a_f
            v = v_f

            #Offset set distance above/below host
            x = a_f * p_a #Impact parameter, multiples of a_p
            y = -10000 #As above, AU
            #z = 0
            r = np.sqrt(x**2 + y**2)
            
            #Components of velocity
            vx = 0 #x-component is always zero
            vy = v_f
    
    
        writeList(body_id,pid,step,time,M,x,y,z,r,vx,vy,vz,v,bodies)
    
    write_args(eccentricity, Mp, a_p, v_p, M_f, Ms, a_f, v_f, ta)
        
    return bodies, n_bodies  



def calc_com(ta,Ms, Mp, p_x, p_y, a, p_vx, p_vy, p_v, eccentricity):
    #Centre of Mass
    com_to_host = (Mp*a/(Ms + Mp))

    com_x = (Mp*p_x/(Ms + Mp))
    com_y = (Mp*p_y/(Ms + Mp))

    com_v = (Mp*p_v)/(Ms + Mp)
            
    com_vx = (Mp*p_vx)/(Ms + Mp)
    com_vy = (Mp*p_vy)/(Ms + Mp)

    return com_to_host,com_x,com_y,com_v,com_vx,com_vy


def change_ref_frame(a,x,y,r,v,vx,vy,com_to_host,com_x,com_y,com_v,com_vx,com_vy):
    
    a -= com_to_host #AU
    
    x -= com_x
    y -= com_y
    z = 0
    r = np.sqrt(x**2+y**2+z**2)
    
    vx -= com_vx
    vy -= com_vy
    vz = 0
    v = np.sqrt(vx**2+vy**2+vz**2)
    
    return a,x,y,z,r,vx,vy,vz,v
    
def write_args(eccentricity, mass, a, v, M_f, Ms, a_f, v_f,ta):
    
    argsFile = open('system_args.dat', "w")
    
    print("Passed arguments of system", file = argsFile)
    
    print("\nMass of star: {}".format(Ms), file = argsFile)
    
    print("\nEccentricity of planet: {}".format(eccentricity), file = argsFile)
    print("\nMass of planet: {}".format(mass), file = argsFile)
    print("\nSemi-major axis of planet: {}".format(a), file = argsFile)
    print("\nVelocity of planet: {}".format(v), file = argsFile)
    
    print("\nMass of Perturber: {}".format(M_f), file = argsFile)
    print("\nImpact Parameter: {}*a".format(a_f), file = argsFile)
    print("\nVelocity of Perturber: {}".format(v_f), file = argsFile)

    ta = ta*(180/np.pi)

    print("\nTrue anomaly: {}".format(ta), file = argsFile)
    


if __name__ == '__main__':
    
    eccentricity = 0
    M_f = 0.00095
    a_f = 2
    v_f = 0.2
    m_p = 1
    p_a = 100
    
    makeIC_file(eccentricity, M_f, a_f, v_f, m_p, p_a)
    
    print('\nDone')
