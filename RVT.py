import numpy as np
import matplotlib.pyplot as plt 
from scipy import special
import time as time_package

start_time = time_package.time()

t = 0
#t_end = 10 * 24 * 3600
t_end = 50000
dt = 1

k = 2
c = 4000
rho = 1000
L = 320000
T_surf = 253.15
T_bot = 273.15

y_m = 0
dx =  0.0005
dy_previous = 0.001
dy_count = 0
time = []
dy_list = []
y_m_list = []

time.append(0)
y_m_list.append(0)


#dy_list.append(0)
dy_list.append(dy_previous / 2)
dy_list.append(dy_previous)

temp_last = [T_surf, 263.15, T_bot]  


def Euler():
    dy = k / (rho * L) * (temp_last[-1] - temp_last[-2]) / dy_previous * dt
    return dy

def Explicit(h):
    temp_new = temp_last[h] + k / (rho * c) * (dt / dy_list[h]**2) * (temp_last[h + 1] - 2 * temp_last[h] + temp_last[h - 1])
    return temp_new

def Explicit_nonuniform(h):
    dy1 = dy_list[h-1]
    dy2 = dy_list[h]
    temp_new = temp_last[h] + k / (rho * c) * dt * (
        2 / (dy1 * (dy1 + dy2)) * (temp_last[h-1] - temp_last[h]) +
        2 / (dy2 * (dy1 + dy2)) * (temp_last[h+1] - temp_last[h])
    ) / 2
    return temp_new



t += dt 
y_m = y_m + dy_previous 
time.append(t)
y_m_list.append(y_m)
dy_count = dy_previous
while t < t_end: 
    dy_previous = Euler() 
    y_m = y_m + dy_previous 
    dy_count += dy_previous
    time.append(t)
    y_m_list.append(y_m)

    if dy_count > dx:
        dy_list.append(dy_count) 
        dy_previous = dy_count
        dy_count = 0
        temp_last.append(T_bot)
    else: 
        dy_list.pop()
        dy_list.append(dy_count)
        dy_previous = dy_count
    
    #print(dy_list)
    #print(temp_last)

    temp_in_this_step = [] 
    temp_in_this_step.append(T_surf) 
    dt_max = 1 * dy_list[-1]**2 * rho * c / k 
    dt = min(dt_max, 1) 
    dt = max(dt, 10**-7)

    for h in range(1, len(dy_list) - 1): 
        #temp_new = Explicit(h) #
        temp_new = Explicit_nonuniform(h)
        temp_in_this_step.append(temp_new) 

    temp_in_this_step.append(T_bot) 
    temp_last = temp_in_this_step 
    
    #print(temp_last)
    # print("dt") 
    #print(dt) 
    # print("cas")
    if int(t) % 1000 == 0:
        print(t) 
    # print("y_m") 
    #print(y_m)
    t += dt 

#print(dy_list)
#print(y_m)
#print(temp_last)

end_time = time_package.time()  # record end

elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time:.2f} seconds")
plt.figure(figsize=(8, 5))
plt.plot(time, y_m_list, label="y_m over time")
plt.xlabel("Time (s)")
plt.ylabel("y_m (m)")
plt.title("Evolution of y_m over Time")
plt.legend()
plt.grid(True)
plt.show()