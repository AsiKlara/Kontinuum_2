import numpy as np
import matplotlib.pyplot as plt 
import scipy


#dT/dt = -k * d^2T/dx^2
#      time, y 

t = 0
t_end = 10 * 24 * 3600
#t_end = 100
tau = 0
h = 0
dt = 1

k = 2
c = 4000
rho = 1000
L = 320000
T_surf = 253.15
T_bot = 273.15

y_m = 0
dx =  0.00005
dy_previous = 0.0001
dy_count = 0
time = []
dy_list = []

#time = list(range(int(t_end)))

dy_list.append(0)
dy_list.append(dy_previous / 2)
dy_list.append(dy_previous)

temp_last = [T_surf, 263.15, T_bot]  


def Euler():
    dy = k / (rho * L) * (temp_last[-1] - temp_last[-2]) / dy_previous * dt
    return dy

def Explicit(h):
    temp_new = temp_last[h] + k / (rho * c) * (dt / dy_previous**2) * (temp_last[h + 1] - 2 * temp_last[h] + temp_last[h - 1])
    return temp_new

# Function for TMDA Algorithm
# https://gist.github.com/vuddameri/75212bfab7d98a9c75861243a9f8f272
def thomas(a, b, c, d):
    """ A is the tridiagnonal coefficient matrix and d is the RHS matrix"""
    N = len(a)
    cp = np.zeros(N,dtype='float64') # store tranformed c or c'
    dp = np.zeros(N,dtype='float64') # store transformed d or d'
    X = np.zeros(N,dtype='float64') # store unknown coefficients
    
    # Perform Forward Sweep
    # Equation 1 indexed as 0 in python
    cp[0] = c[0]/b[0]  
    dp[0] = d[0]/b[0]
    # Equation 2, ..., N (indexed 1 - N-1 in Python)
    for i in np.arange(1,(N),1):
        dnum = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/dnum
        dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
    # Perform Back Substitution
    X[(N-1)] = dp[N-1]  # Obtain last xn 

    for i in np.arange((N-2),-1,-1):  # use x[i+1] to obtain x[i]
        X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
    return X

def Crank_Nicolson():
    a = []
    b = []
    c = []
    d = []
    b.append(1.0)
    c.append(0.0)
    d.append(T_surf)

    print(temp_last)
    for h in range(len(dy_list) - 1):
        print(h)
        mu = (dt / dy_list[h + 1]**2)
        a.append(1/2 * mu) 
        b.append(1 + mu)
        c.append(1 / 2 * mu)
        d.append(1 + 1/2 * mu * (temp_last[h + 1] - 2 * temp_last[h] + temp_last[h - 1]))

    a.append(0.0)
    b.append(1.0)
    d.append(T_bot)

    # Solve tridiagonal system
    temp_in_this_step = thomas(a[1:], b, c[:-1], d)
    #temp_in_this_step = thomas(a, b, c, d)
    return temp_in_this_step
    #temp_new = temp_last[h] + (dt / dy_list[h + 1]**2) * (temp_last[h + 1] - 2 * temp_last[h] + temp_last[h - 1])


t += dt 
y_m = y_m + dy_previous 
while t < t_end: 
    dy_previous = Euler() 
    y_m = y_m + dy_previous 
    dy_count += dy_previous
    if dy_count > dx:
        dy_list.append(dy_count) 
        dy_previous = dy_count
        dy_count = 0
    else: 
        dy_list.pop()
        dy_list.append(dy_count)
        dy_previous = dy_count
    temp_in_this_step = [] 
    temp_in_this_step.append(T_surf) 
    dt_max = 0.5 * dy_previous**2 * rho * c / k 
    dt = min(dt_max, 1) 
    for h in range(len(dy_list) - 1): #
        temp_new = Explicit(h) #
        temp_in_this_step.append(temp_new) 
        #temp_in_this_step = Crank_Nicolson() 
        # #print(h) 
        temp_in_this_step.append(T_bot) 
        temp_last = temp_in_this_step 
        print("dt") 
        print(dt) 
        print("cas")
        print(t) 
        print("y_m") 
        print(y_m)
        t += dt 
        #print(y_m)