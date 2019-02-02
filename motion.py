import matplotlib.pyplot as plt
import numpy as np


m1, m2 = 15., 15.
l1, l2 = 5., 5.
g = 9.8

def DiffEq(y, t):
    b = (2.*m1+m2-m2*np.cos(2.*y[0,0]-2.*y[1,0]))
    eq = np.zeros((2,2))
    eq[0,0] = y[0,1]
    eq[1,0] = y[1,1]
    eq[0,1] = (-g*(2.*m1+m2)*np.sin(y[0,0]) - 
		m2*g*np.sin(y[0,0]-2*y[1,0]) -
               	2.*np.sin(y[0,0]-y[1,0])*m2*(y[1,1]**2.*l2 + 
                y[0,1]**2.*l1*np.cos(y[0,0]-y[1,0])))/(l1*b)
    
    eq[1,1] = (2.*np.sin(y[0,0]-y[1,0])*(y[0,1]**2.*l1*(m1+m2)+g*(m1+m2)*np.cos(y[0,0])+
               y[1,1]**2.*l2*m2*np.cos(y[0,0]-y[1,0])))/(l2*b)
    return eq


def start():  
    theta1, theta2 = np.pi/2., np.pi/2.
    Y = np.array([[theta1,0.0],[theta2,0.0]])	#intital values: theta1, vel_theta1, theta2, vel_theta2
    t, h,t_tot = 0.0, 0.01, 50. 
    time = []
    ang1, ang2 = [], []
    x_ball1,y_ball1, ball1 = [], [], []
    x_ball2,y_ball2, ball2 = [], [], []

    while(t < t_tot): 
        Y = RK4(DiffEq, Y, t, h) 
        time.append(t)
        ang1.append(np.rad2deg(Y[0,0])%360.)
        ang2.append(np.rad2deg(Y[1,0])%360.)
        x1, y1 = l1*np.sin(Y[0,0]), -l1*np.cos(Y[0,0])
        x2, y2 = x1 + l2*np.sin(Y[1,0]), y1 - l1*np.cos(Y[1,0])
        x_ball1.append(x1)
        y_ball1.append(y1)
        x_ball2.append(x2)
        y_ball2.append(y2)
        print("angle1:", np.rad2deg(Y[0,0])%360.)
        print("angle2:", np.rad2deg(Y[1,0])%360.)
        t = t + h

    #the x and y of ball 1 and ball 2 
    plt.figure()
    plt.plot(x_ball1, y_ball1)
    plt.plot(x_ball2, y_ball2)
    plt.show()
    #x and y of ball 1 with respect of time
    plt.figure()
    plt.plot(time, x_ball1)
    plt.plot(time, y_ball1)
    plt.show()
    #x and y of ball 2 with respect of time
    plt.figure()
    plt.plot(time, x_ball2)
    plt.plot(time, y_ball2)
    plt.show()
    #angle  of ball 1 and 2 with respect to time 
    plt.figure()
    plt.plot(time, ang1)
    plt.plot(time, ang2)
    plt.show()

def RK4(diffeq, y0, t, h):
    """ RK4 method for ODEs:
        Given y0 at t, returns y1 at t+h 
        For a system of ODEs, diffeq must return a NumPy array"""
    k1 = h*diffeq(y0, t)                    # dy/dt at t
    k2 = h*diffeq(y0+0.5*k1, t + h/2.)      # dy/dt at t+h/2
    k3 = h*diffeq(y0+0.5*k2, t + h/2.)      # dy/dt at t+h/2
    k4 = h*diffeq(y0+k3, t + h)             # dy/dt at t+h
    return y0 + (k1+k4)/6.0 + (k2+k3)/3.0    

start()
