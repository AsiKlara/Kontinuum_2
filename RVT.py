import numpy as np
import matplotlib.pyplot as plt 
from scipy import special
import time as time_package

start_time = time_package.time()

t = 0
t_end = 10 * 24 * 3600
#t_end = 100000
dt = 1

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
dx =  0.001
dy_previous = 0.002
# dx = 0.0005
# dy_previous = 0.001
dy_count = 0
time = []
dy_list = []
y_m_list = []

# for plotting
temp_profiles = []  
temp_profiles_stefan = []
y_profiles = []     

time.append(0)
y_m_list.append(0)


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


def lambda_1(lmbda):
    return lmbda * np.exp(lmbda**2) * special.erf(lmbda) - c * (T_surf - T_bot) / (L  * np.sqrt(np.pi))

def calculate_lambda(lmbda):
    lmbda2 = lambda_1(lmbda)
    while abs(lmbda - lmbda2) / abs(lmbda2) < 0.00001:
        lmbda = lmbda2
        lmbda2 = lambda_1(lmbda2)
    lmbda = lmbda2
    print(lmbda)
    return lmbda

def stefan_y_m(time):
    return 2 * lmbda2 * np.sqrt(kappa * time)

def stefan_temp(lmbda2):
    temp_stefan = []
    y = dy_list[0]
    for i in range(1, len(dy_list)):
        temp_stefan.append(special.erf(y / (2 * np.sqrt(kappa * t))) / special.erf(lmbda2))
        y += dy_list[i]
    return temp_stefan

lmbda2 = calculate_lambda(lmbda)

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

    temp_in_this_step = [] 
    temp_in_this_step.append(T_surf) 
    dt_max = 1 * dy_list[-1]**2 * rho * c / k 
    #dt = min(dt_max, 1) 
    dt = max(dt_max, 10**-7)

    for h in range(1, len(dy_list) - 1): 
        temp_new = Explicit_nonuniform(h)
        temp_in_this_step.append(temp_new) 

    temp_in_this_step.append(T_bot) 
    temp_last = temp_in_this_step 

    if int(t) % 100 == 0:
        print(t) 
        #print(temp_in_this_step)
        y_positions = np.cumsum(dy_list)
        y_profiles.append(y_positions)
        temp_profiles.append(temp_in_this_step)
        temp_stefan = stefan_temp(lmbda2)
        temp_profiles_stefan.append(temp_stefan)

    t += dt 

y_m_stefan = []
for i in range(len(time)):
    y_m_stefan.append(stefan_y_m(time[i]))

end_time = time_package.time()  # record end

elapsed_time = end_time - start_time
print(f"Simulation took {elapsed_time:.2f} seconds")
plt.figure(figsize=(8, 5))
plt.plot(time, y_m_list, label="Numerical y_m")
plt.plot(time, y_m_stefan, label="Analytical y_m")
plt.xlabel("Time (s)")
plt.ylabel("y_m (m)")
plt.title("Evolution of y_m over Time")
plt.legend()
plt.grid(True)
plt.show()

# this plotting is AI generated, because Im really lazy
# --- Build 2D color map of temperature field over time ---
max_depth = max(y[-1] for y in y_profiles)
ny = 300  # more resolution in depth
y_grid = np.linspace(0, max_depth, ny)

# Create uniform time grid corresponding to saved profiles
time_array = np.linspace(0, t_end, len(temp_profiles))

# Interpolate numerical freezing front to this time grid
y_m_interp = np.interp(time_array, time, y_m_list)

# Initialize temperature matrix
temp_matrix = np.full((ny, len(temp_profiles)), np.nan)

# --- Robust interpolation for nonuniform y_profiles ---
for i, (y, T) in enumerate(zip(y_profiles, temp_profiles)):
    # Ensure y and T are arrays and sorted (sometimes dy_list changes)
    y = np.array(y)
    T = np.array(T)
    sort_idx = np.argsort(y)
    y = y[sort_idx]
    T = T[sort_idx]

    # Clip to prevent extrapolation errors
    y_clipped = np.clip(y_grid, y[0], y[-1])

    # Interpolate safely
    temp_matrix[:, i] = np.interp(y_clipped, y, T)

# --- Mask region above the moving front (unfrozen) ---
Y, Tmesh = np.meshgrid(time_array, y_grid)
mask = np.zeros_like(temp_matrix, dtype=bool)
for i, ym in enumerate(y_m_interp):
    mask[:, i] = y_grid >= ym  # above front = unfrozen

temp_masked = np.ma.array(temp_matrix, mask=mask)

# --- Make colormap with white for masked values ---
cmap = plt.cm.coolwarm.copy()
cmap.set_bad(color='white')

# --- Plot ---
plt.figure(figsize=(9, 5), facecolor='white')
plt.pcolormesh(time_array, y_grid, temp_masked, shading='auto', cmap=cmap)
plt.colorbar(label="Temperature (K)")
plt.xlabel("Time (s)")
plt.ylabel("Depth (m)")
plt.title("Temperature evolution")
plt.gca().invert_yaxis()
plt.grid(False)
plt.tight_layout()
plt.show()


# for i in range(len(temp_profiles)):
#     plt.plot(temp_profiles[i], y_profiles[i], label=f"t = {int(i*1000)} s")

# plt.gca().invert_yaxis()  # surface at top (y=0)
# plt.xlabel("Temperature (K)")
# plt.ylabel("Depth (m)")
# plt.title("Temperature profiles during freezing")
# plt.legend()
# plt.grid(True)
# plt.show()