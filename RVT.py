import numpy as np
import matplotlib.pyplot as plt 
from scipy import special
import time as time_package

start_time = time_package.time()

t = 0
t_end = 3 * 24 * 3600
#t_end = 100000

# odhad lambda
lmbda = 0.5
k = 2
c = 4000
rho = 1000
kappa = k / (rho  * c)
L = 320000
T_surf = 253.15
T_bot = 273.15

y_m = 0
dx =  0.002
dy_previous = 0.004
# dx = 0.0005
# dy_previous = 0.001
# dx = 0.005
# dy_previous = 0.01
dy_count = 0

dt = 0.5 * dx**2 * rho * c / k 


time = []
time_plot_list = []
dy_list = []
y_m_list = []

# for plotting
temp_profiles = []  
temp_profiles_stefan = []
y_profiles = []  
y_profiles_stefan = []   
y_stefan = []

time.append(0)
y_m_list.append(0)
y_stefan.append(0)


dy_list.append(dy_previous / 2)
dy_list.append(dy_previous / 2)

temp_last = [T_surf, 263.15, T_bot]  
temp_1 = [T_surf, 263.15, T_bot]  
dy_1 = []
dy_1.append(0)
dy_1.append(dy_previous / 2)
dy_1.append(dy_previous)

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
    ) 
    return temp_new


def lambda_1(lmbda):
    return lmbda * np.exp(lmbda**2) * special.erf(lmbda) - c * (T_surf - T_bot) / (L  * np.sqrt(np.pi))


def calculate_lambda(lmbda):
    lmbda2 = lambda_1(lmbda)
    while abs(lmbda - lmbda2) / abs(lmbda2) < 0.00001:
        lmbda = lmbda2
        lmbda2 = lambda_1(lmbda2)
    lmbda = lmbda2
    return lmbda


def stefan_y_m(time):
    return 2 * lmbda2 * np.sqrt(kappa * time)


def stefan_temp(lmbda2):
    temp_stefan = []
    y = y_stefan[0]
    for i in range(1, len(y_stefan)):
        temp_stefan.append(special.erf(y / (2 * np.sqrt(kappa * t))) / special.erf(lmbda2))
        y = y_stefan[i]
    return temp_stefan


lmbda2 = calculate_lambda(lmbda)

t += dt 
y_m = y_m + dy_previous 
time.append(t)
y_m_list.append(y_m)
dy_count = 0
appended = True
time_save_previous = 0
step = 0

while t < t_end: 
    dy_previous = Euler() 
    y_m = y_m + dy_previous 
    dy_count += dy_previous
    time.append(t)
    y_m_list.append(y_m)

    if dy_count >= dx:
        dy_list[-1] = dx
        temp_last.append(T_bot)  
        dy_previous = dx
        dy_count -= dx
        appended = True
    else:
        if appended:
            dy_list.append(dy_count)
            temp_last.append(T_bot) 
            appended = False
        else:
            dy_list[-1] = dy_count
        dy_previous = dy_count

    temp_in_this_step = [] 
    temp_in_this_step.append(T_surf) 
    #dt_max = 0.5 * dy_list[-1]**2 * rho * c / k 
    dt_max = 0.1 * dx**2 * rho * c / k 
    dt = max(dt_max, 10**-7)

    for h in range(1, len(dy_list) - 1): 
        temp_new = Explicit_nonuniform(h)
        temp_in_this_step.append(temp_new) 

    temp_in_this_step.append(T_bot) 
    temp_last = temp_in_this_step 

    if int(t) % 1000 == 0:
        if time_save_previous != int(t):
            # print(t) 
            time_plot_list.append(t / 3600)
            y_positions = np.cumsum(dy_list)
            y_profiles.append(y_positions)
            y_m_stef = float(stefan_y_m(t))
            temp_profiles.append(temp_in_this_step)
            n_points = len(dy_list)
            
            y_stefan_profile = np.linspace(0, y_m_stef, n_points)

            temp_stefan_profile = 253.15 + 20 * (
                special.erf(y_stefan_profile / (2 * np.sqrt(kappa * t))) / special.erf(lmbda2)
            )

            y_profiles_stefan.append(y_stefan_profile)
            temp_profiles_stefan.append(temp_stefan_profile)
        time_save_previous = int(t)

    t += dt 

temp_stefan_last = 253.15 + 20 *  (special.erf(y_stefan_profile / (2 * np.sqrt(kappa * t))) / special.erf(lmbda2))

y_m_stefan = []
for i in range(len(time)):
    y_m_stefan.append(stefan_y_m(time[i]))

end_time = time_package.time() 

elapsed_time = end_time - start_time

for i in range(len(time)):
    time[i] = time[i] / 3600


print(f"Simulation took {elapsed_time:.2f} seconds")
plt.figure(figsize=(8, 5))
plt.plot(time, y_m_list, label="Numerical y_m")
plt.plot(time, y_m_stefan, label="Analytical y_m")
plt.xlabel("Time (h)")
plt.ylabel("y_m (m)")
plt.title("Evolution of y_m over Time")
plt.legend()
plt.grid(True)
plt.savefig("y_m5.png")

with open("y.txt", "w") as f:
    for i in range(len(time)):
        f.write(f"{time[i]:.5f}    {y_m_list[i]:.5f}    {y_m_stefan[i]:.5f}\n")


y_numerical = []
y_numerical.append(0)
y_sum = 0
for i in range(len(dy_list) - 1):
    y_sum += dy_list[i]
    y_numerical.append(y_sum)


plt.figure(figsize=(8, 5))
plt.plot()
plt.plot(temp_1, dy_1, 'o', label="Initial", linestyle='none')
plt.plot(temp_last, y_numerical, '^', label="Numerical", linestyle='none')
plt.plot(temp_stefan_last, y_stefan_profile, '*', label="Analytical", linestyle='none')
plt.xlabel("Temperature (K)")
plt.ylabel("y_m (m)")
plt.gca().invert_yaxis()
plt.title("Evolution of y_m over Time")
plt.legend()
plt.grid(True)
plt.savefig("temperatre_in_ice5.png")

with open("initial.txt", "w") as f:
    for i in range(len(dy_1)):
        f.write(f"{dy_1[i]:.5f}    {temp_1[i]:.5f}\n")

with open("numerical.txt", "w") as f:
    for i in range(len(y_numerical)):
        f.write(f"{y_numerical[i]:.5f}    {temp_last[i]:.5f}\n")

with open("analytical.txt", "w") as f:
    for i in range(len(y_stefan_profile)):
        f.write(f"{y_stefan_profile[i]:.5f}    {temp_stefan_last[i]:.5f}\n")

# this plotting is AI generated, because Im really lazy
# --- Scatter plot of raw temperature profiles over time ---
def plot_temp(y_profiles, temp_profiles):

    # Collect all points
    t_points = []
    y_points = []
    T_points = []

    for t_val, y, T in zip(time_plot_list, y_profiles, temp_profiles):
        y = np.array(y)
        T = np.array(T)
        # repeat current time for all depths at this snapshot
        t_points.extend([t_val] * len(y))
        y_points.extend(y)
        T_points.extend(T)

    # Convert to numpy arrays
    t_points = np.array(t_points)
    y_points = np.array(y_points)
    T_points = np.array(T_points)

    # --- Scatter plot ---
    sc = plt.scatter(
        t_points,
        y_points,
        c=T_points,
        cmap='coolwarm',
        s=10,       # dot size
        edgecolors='none'
    )

    # Add colorbar and labels
    plt.colorbar(sc, label="Temperature (K)")
    plt.xlabel("Time (h)")
    plt.ylabel("Depth (m)")
    plt.gca().invert_yaxis()  # deeper = down
    plt.tight_layout()


plt.figure(figsize=(9, 5), facecolor='white')
plt.title("Temperature evolution Numerical")
plot_temp(y_profiles, temp_profiles)
plt.savefig("temperature5.png")


plt.figure(figsize=(9, 5), facecolor='white')
plt.title("Temperature evolution Analytical")
plot_temp(y_profiles_stefan, temp_profiles_stefan)
plt.savefig("temperature_stefan5.png")