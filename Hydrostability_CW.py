import numpy as np
import matplotlib.pyplot as plt

def get_PE_GZ_case1(theta, rho):
# For the case where the submerged part is a trapezium
# rho is the RELATIVE density of the solid to the water
    
    t = np.tan(theta)
    c = np.cos(theta)
    s = np.sin(theta)

    i = np.array([1,0,0])
    j = np.array([0,1,0])
    k = np.array([0,0,1])

    beta = rho - t/2
    Arect = beta
    Atri = t/2
    Awater = Arect + Atri

    y_bar = (Arect*(beta/2) + Atri*(beta+t/3)) / (Awater)
    x_bar = (1/6) * ((y_bar - beta/2)/(beta/2 + t/3))

    # Define centre of bouyancy (CoG of displaced water)
    CoB = x_bar*i + y_bar*j

    # Define centre of gravity (CoG of solid)
    CoG = 0.5*j

    delta = CoG - CoB

    # Define unit vectors v,u wrt the water surface
    v = -s*i + c*j
    u = c*i + s*j

    PE = rho*9.807*(np.dot(delta,v))
    GZ = -np.dot(delta,u)
    return([PE,GZ,CoB])


def get_PE_GZ_case2(theta, rho):
# For the case where the submerged part is a triangle
# rho is the RELATIVE density of the solid to the water
    
    t = np.tan(theta)
    c = np.cos(theta)
    s = np.sin(theta)

    i = np.array([1,0,0])
    j = np.array([0,1,0])
    k = np.array([0,0,1])

    base = np.sqrt((2*rho)/(t))
    height = base * t

    y_bar = height/3
    x_bar = 0.5 - base/3

    # Define centre of bouyancy (CoG of displaced water)
    CoB = x_bar*i + y_bar*j

    # Define centre of gravity (CoG of solid)
    CoG = 0.5*j

    delta = CoG - CoB

    # Define unit vectors v,u wrt the water surface
    v = -s*i + c*j
    u = c*i + s*j

    PE = rho*9.807*(np.dot(delta,v))
    GZ = -np.dot(delta,u)
    return([PE,GZ,CoB])

def get_PE_GZ_case3(theta, rho):
# For the case where the unsubmerged part is a triangle
# rho is the RELATIVE density of the solid to the water
    
    t = np.tan(theta)
    c = np.cos(theta)
    s = np.sin(theta)

    i = np.array([1,0,0])
    j = np.array([0,1,0])
    k = np.array([0,0,1])

    base = np.sqrt((2*(1-rho))/(t))
    height = base * t

    y_bar = (0.5 - (1-(height/3))*(1-rho))/rho
    x_bar = -(((base/3)-0.5)*(1-rho))/rho

    # Define centre of bouyancy (CoG of displaced water)
    CoB = x_bar*i + y_bar*j

    # Define centre of gravity (CoG of solid)
    CoG = 0.5*j

    delta = CoG - CoB

    # Define unit vectors v,u wrt the water surface
    v = -s*i + c*j
    u = c*i + s*j

    PE = rho*9.807*(np.dot(delta,v))
    GZ = -np.dot(delta,u)
    return([PE,GZ,CoB])

def get_PE(theta, rho):
    if (2*rho < np.tan(theta)) and (rho < 0.5):
        PE = get_PE_GZ_case2(theta, rho)[0]
        return(PE)
    if (2 - 2*rho < np.tan(theta)) and (rho > 0.5):
        PE = get_PE_GZ_case3(theta, rho)[0]
        return(PE)
    else:
        PE = get_PE_GZ_case1(theta,rho)[0]
        return(PE)
    

def plot_PE_against_angle_for_given_rho(rho):
    heel_angles = np.linspace(0, np.pi*45/180, 600)
    PEs = np.zeros(600)
    for index in range(len(heel_angles)):
        PEs[index] = get_PE(heel_angles[index], rho)

    plt.plot(heel_angles*(180/np.pi), PEs)
    min_PE_angle = heel_angles[np.argmin(PEs)]
    plt.axvline(x = min_PE_angle*(180/np.pi), color = 'orange')
    print(min_PE_angle*(180/np.pi))
    plt.grid()
    plt.title('Potential energy against heel angle for constant realtive density:\n(orange line shows stable angle)')
    plt.xlabel('Heel Angle, degrees')
    plt.ylabel('Potential Energy')
    plt.show()

def plot_relative_density_against_stable_angle():
    rho_list = np.linspace(0.1,0.99,99)
    min_PE_angles = np.zeros(99)
    for n in range(len(rho_list)):

        heel_angles = np.linspace(0, (45*np.pi)/180, 600)
        PEs = np.zeros(600)
        for m in range(len(heel_angles)):
            PEs[m] = get_PE(heel_angles[m], rho_list[n])

        min_PE_angles[n] = heel_angles[np.argmin(PEs)]
    plt.plot(min_PE_angles*(180/np.pi),rho_list)
    plt.plot(-min_PE_angles*(180/np.pi),rho_list)
    plt.axis([-50, 50, 0, 1])
    plt.grid()
    plt.title('Bifurcation diagram for stable equilibrium')
    plt.xlabel('Stable heel angle, degrees')
    plt.ylabel('Relative density')
    plt.show()

def get_GZ(theta, rho):
    if (2*rho < np.tan(theta)) and (rho < 0.5):
        GZ = get_PE_GZ_case2(theta, rho)[1]
        return(GZ)
    if (2 - 2*rho < np.tan(theta)) and (rho > 0.5):
        GZ = get_PE_GZ_case3(theta, rho)[1]
        return(GZ)
    else:
        GZ = get_PE_GZ_case1(theta,rho)[1]
        return(GZ)

def plot_GZ_against_heel_angle_for_given_rho(rho):
    heel_angles = np.linspace(0, 45*np.pi/180, 600)
    GZs = np.zeros(600)
    for index in range(len(heel_angles)):
        GZs[index] = get_GZ(heel_angles[index], rho)

    plt.plot(heel_angles*(180/np.pi), GZs)
    plt.grid()
    plt.title('GZ Method for constant relative density')
    plt.xlabel('Heel angle, degrees')
    plt.ylabel('GZ (+ve GZ acts to reduce heel angle)')
    plt.show()

def get_CoB(theta, rho):
    if (2*rho < np.tan(theta)) and (rho < 0.5):
        CoB = get_PE_GZ_case2(theta, rho)[2]
        return(CoB)
    if (2 - 2*rho < np.tan(theta)) and (rho > 0.5):
        CoB = get_PE_GZ_case3(theta, rho)[2]
        return(CoB)
    else:
        CoB = get_PE_GZ_case1(theta,rho)[2]
        return(CoB)

def catastrophe_theory_plot(rho):
    heel_angles = np.linspace(0, (45*np.pi)/180, 600)
    CoBs = np.zeros([600,3])
    for i in range(len(heel_angles)):
        CoBs[i] = get_CoB(heel_angles[i], rho)
    CoBs_x = np.zeros(600)
    CoBs_y = np.zeros(600)
    for i in range(len(CoBs)):
        CoBs_x[i] = CoBs[i,0]
        CoBs_y[i] = CoBs[i,1]
    
    fig, ax = plt.subplots()
    plt.plot(CoBs_x,CoBs_y, color = 'blue', label = 'Bouyancy locus')
    plt.plot(-CoBs_x,CoBs_y, color = 'blue')

    verticals = np.zeros([600,2])
    i = np.array([1,0])
    j = np.array([0,1])
    for index in range(len(heel_angles)):
        theta = heel_angles[index]
        c = np.cos(theta)
        s = np.sin(theta)
        v = -s*i + c*j # <-- what we're after
        verticals[index] = v
    
    for i in range(len(CoBs)//35):
        i = i*35
        x_for_line = np.array([CoBs[i,0],CoBs[i,0]+verticals[i,0]])
        y_for_line = np.array([CoBs[i,1],CoBs[i,1]+verticals[i,1]])
        plt.plot(x_for_line,y_for_line,color = 'orange', linewidth = 0.7)
        plt.plot(-x_for_line,y_for_line,color = 'orange', linewidth = 0.7)
    
    plt.plot([0],[0.5],'o', markersize = 3, color = 'black', label = 'Centre of Gravity')
    plt.legend()
    plt.title('Construction of the evolute for relative density = 0.4')
    plt.xlabel('x')
    plt.ylabel('y')
    ax.set_aspect('equal')
    plt.show()

# plot_relative_density_against_stable_angle()
plot_PE_against_angle_for_given_rho(0.)

plot_GZ_against_heel_angle_for_given_rho(0.4)

catastrophe_theory_plot(0.4)

""" UNCOMMENT SOME/ALL OF THE ABOVE TO GET PLOTS """