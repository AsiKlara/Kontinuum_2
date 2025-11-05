import numpy as np
import matplotlib.pyplot as plt 
from scipy import special
import time as time_package

t = 0
t_end = 10 * 24 * 3600
t_end = 700000
dt = 1

k = 2
c = 4000
rho = 1000
kappa = k / (rho  * c)
L = 320000
T_surf = 253.15
T_bot = 273.15

y_m = 10
time = []
time_plot_list = []
dy_list = []
y_m_list = []

# for plotting
temp_profiles = []    

time.append(0)

dy_list = np.ones(100) * 0.01

temp_last = np.ones(100) * 273.15

temp_profiles.append(temp_last)
time_plot_list.append(0)

def Explicit_nonuniform(h):
    dy1 = dy_list[h-1]
    dy2 = dy_list[h]
    temp_new = temp_last[h] + k / (rho * c) * dt * (
        2 / (dy1 * (dy1 + dy2)) * (temp_last[h-1] - temp_last[h]) +
        2 / (dy2 * (dy1 + dy2)) * (temp_last[h+1] - temp_last[h])
    ) 
    return temp_new

t += dt 
time.append(t)
t_print = 0
while t < t_end: 
    #print(temp_last)
    time.append(t)

    temp_in_this_step = [] 
    temp_in_this_step.append(T_surf) 
    dt_max = 0.5 * dy_list[-1]**2 * rho * c / k 
    dt = max(dt_max, 10**-7)
    dt = min(dt_max, 100)

    for h in range(1, len(dy_list) - 1): 
        temp_new = Explicit_nonuniform(h)
        temp_in_this_step.append(float(temp_new)) 

    temp_in_this_step.append(T_bot) 
    temp_last = temp_in_this_step 


    if int(t_print) / 50000 >= 1:
        time_plot_list.append(int(t) / (3600))
        temp_profiles.append(temp_last)
        t_print = 0

    t_print += dt
    t += dt 

#print(len(temp_last))
# print(len(temp_profiles[0]))

# --- Build y positions from dy_list (cumulative depth from surface) ---
y_list = np.cumsum(dy_list)
# y_list = np.insert(y_list, 0, 0)
# print(len(y_list))

# --- Plot temperature profiles at saved times ---
plt.figure(figsize=(8, 5), facecolor='white')

for i in range(len(temp_profiles)):
    #print(temp_profiles[i])
    plt.plot(temp_profiles[i], y_list, label=f"t = {int(time_plot_list[i])} hours")

plt.gca().invert_yaxis()  
plt.xlabel("Temperature (K)")
plt.ylabel("Depth (m)")
plt.title("Temperature profiles at different times")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("rvt_test.pdf")