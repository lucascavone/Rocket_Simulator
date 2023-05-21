import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

g = 9.8101 # m/s^2
air_density = 1.225 # kg/m^3

# When called, prompts user for float value larger than 0
def get_inputs(input_name):
    while True:

        while True:
            input_value = input(f"Enter value for {input_name}:\n >> ")
            try:
                input_value = float(input_value)
                break
            except ValueError:
                print("Enter an appropriate value.")

        if input_value > 0:
            break
        else:
            print(f"{input_name} cannot be negative or equal to 0.")

    return input_value

def get_initial_height(input_name):
    while True:
        while True:
            input_value = input(f"Enter value for {input_name}:\n >> ")
            try:
                input_value = float(input_value)
                break
            except ValueError:
                print("Enter an appropriate value.")

        if input_value >= 0:
            break
        else:
            print(f"{input_name} cannot be negative.")
    return input_value

def get_radius(input_name):
    while True:
        while True:
            input_value = input(f"Enter value for {input_name}:\n >> ")
            try:
                input_value = float(input_value)
                break
            except ValueError:
                print("Enter an appropriate value.")

        if input_value > 0:
            area = input_value ** 2 * math.pi
            break
        else:
            print(f"{input_name} cannot be negative or equal to 0.")
    return area

# y vs time with basic kinematics
def f(yi, vi, a, t):
    yf = yi + vi * t - 0.5 * a * t ** 2
    return yf

# y vs time with air resistance formulas
def h(yi, vi, a, t, vterm):
    yf = yi + (vterm / a) * (vi + vterm) * (1 - math.e ** ((-a * t) / vterm)) - vterm * t
    return yf

# y vs time while engine is firing
def y(mi, a, u, alpha, t, tlift):
    alpha_u = alpha * u
    m = mi - alpha * t
    mlift = mi - alpha * tlift
    yf = ((u ** 2 / a) * (((m * a) / alpha_u) * np.log((m * a) / alpha_u) - (m * a) / alpha_u) -
        0.5 * a * t ** 2 + (mi * a * t) / alpha - u * t) - ((u ** 2 / a) * (((mlift * a) / alpha_u) *
        np.log((mlift * a) / alpha_u) - (mlift * a) / alpha_u) - 0.5 * a * tlift ** 2 +
        (mi * a * tlift) / alpha - u * tlift)
    yf = np.where(yf < 0, 0, yf)
    return yf

while True:
    plot_projectiles_only = input(f"Plot only projectile motion (with & without air resistance)?\n Yes (Y) / No (N) >> ")
    # Only plots projectile motion without rocket motion equations
    if plot_projectiles_only.lower() == "y":
        y_initial = get_initial_height("Initial position (m)")
        v_initial = get_inputs("Initial speed (m/s)")
        projectile_mass = get_radius("Projectile mass (kg)")
        area = get_radius("Projectile radius (m)")
        drag_c = get_inputs("Drag coefficient (~0.2 to ~2, rockets ~ 0.75)")

        # ================BASIC KINEMATICS================
        kinematics_max_height = (v_initial ** 2) / (2 * g) + y_initial

        # Flight time calculation for kinematic trajectory w/o air resistance
        if y_initial == 0:
            flight_time = (2 * v_initial) / g
        elif y_initial > 0:
            flight_time = (v_initial + math.sqrt((v_initial ** 2) + 2 * g * y_initial)) / g

        # ================KINEMATICS WITH AIR RESISTANCE================
        v_term = math.sqrt((2 * projectile_mass * g) / (air_density * area * drag_c))

        # Solving for h(xi, vi, a, t, vterm) = 0 to find flight time
        find_time_air_resis = lambda t: h(y_initial, v_initial, g, t, v_term)
        result_air_resis = fsolve(find_time_air_resis, x0=flight_time)
        root_air_resis = result_air_resis[0]

        # Max height of trajectory calculation
        air_resistance_max_height_time = (-v_term / g) * math.log(v_term / (v_initial + v_term), math.e)
        air_resistance_max_height = h(y_initial, v_initial, g, air_resistance_max_height_time, v_term)

        # ================STATS================
        print(f"\n=============BASIC KINEMATICS============")
        print(f"Max height : {round(kinematics_max_height, 2)}m")
        print(f"Flight time : {round(flight_time, 2)}s")
        print(f"\n===========WITH AIR RESISTANCE===========")
        print(f"Max height : {round(air_resistance_max_height, 2)}m")
        print(f"Flight time : {np.round(root_air_resis, 2)}s")
        print(f"Terminal velocity : {round(v_term, 2)} m/s")

        # ================GRAPHING================
        # Define arrays of time values to plot
        kin_time = np.linspace(0, flight_time, 100) if flight_time else np.linspace(0, 5, 100)
        air_time = np.linspace(0, root_air_resis, 100)

        plt.plot(kin_time, f(y_initial, v_initial, g, kin_time), color='red', label='Kinematic Trajectory')
        plt.plot(air_time, h(y_initial, v_initial, g, air_time, v_term), color='green', label='Kinematics w/ Air Resistance')
        plt.axhline(y=0, color='black')
        plt.legend()
        plt.plot()
        plt.xlabel("Time")
        plt.ylabel("Height")
        plt.show()
        sys.exit()
    elif plot_projectiles_only.lower() == "n":
        break
    else:
        print("Enter an appropriate value.")

# Initialize initial height
y_initial = get_initial_height("Initial position (m)")
# Initialize initial velocity
v_initial = get_inputs("Initial speed (m/s)")
# Initialize mass of rocket
rocket_mass = get_inputs("Total mass (kg)")
# Initialize mass of fuel
fuel_mass = get_inputs("Fuel mass (kg)")
# Initialize cross-sectional area
area = get_radius("Rocket radius (m)")
# Initialize drag coefficient
drag_c = get_inputs("Drag coefficient (~0.75 for rockets)")
# Initialize exhaust speed (u)
exhaust_speed = get_inputs("Exhaust speed (m/s)")
# Initialize exhaust rate (alpha)
exhaust_rate = get_inputs("Exhaust rate (kg/s)")

# ================BASIC KINEMATICS================
kinematics_max_height = (v_initial ** 2) / (2 * g) + y_initial

# Flight time calculation for kinematic trajectory w/o air resistance
if y_initial == 0:
    flight_time = (2 * v_initial) / g
elif y_initial > 0:
    flight_time = (v_initial + math.sqrt((v_initial ** 2) + 2 * g * y_initial)) / g

# ================KINEMATICS WITH AIR RESISTANCE================
v_term = math.sqrt((2 * rocket_mass * g) / (air_density * area * drag_c))

# Solving for h(xi, vi, a, t, vterm) = 0 to find flight time
find_time_air_resis = lambda t: h(y_initial, v_initial, g, t, v_term)
result_air_resis = fsolve(find_time_air_resis, x0=flight_time)
root_air_resis = result_air_resis[0]

# Max height of trajectory calculation
air_resistance_max_height_time = (-v_term / g) * math.log(v_term / (v_initial + v_term), math.e)
air_resistance_max_height = h(y_initial, v_initial, g, air_resistance_max_height_time, v_term)

# ================ROCKET MOTION================
m_dry = rocket_mass - fuel_mass
t_burn = fuel_mass / exhaust_rate
t_liftoff = (rocket_mass / exhaust_rate) - (exhaust_speed / g)
m_liftoff = (exhaust_speed * exhaust_rate) / g

# Velocity at engine cutoff calculation
v_burn = -exhaust_speed * math.log(math.fabs((g * m_dry) / (exhaust_speed * exhaust_rate)), math.e) - g * (t_burn - t_liftoff)

# Height at engine cutoff calculation
y_burn = y(rocket_mass, g, exhaust_speed, exhaust_rate, t_burn, t_liftoff)

# Trajectory calculation after engine cutoff w/o air resistance
rocket_kin_max_height = (v_burn ** 2) / (2 * g) + y_initial + y_burn
rocket_kin_flight_time = (v_burn + math.sqrt((v_burn ** 2) + 2 * g * (y_initial + y_burn))) / g

# Trajectory calculation after engine cutoff w/ air resistance
rocket_vterm = math.sqrt((2 * m_dry * g) / (air_density * area * drag_c))
rocket_air_max_height_time = (-rocket_vterm / g) * math.log(rocket_vterm / (v_burn + rocket_vterm), math.e)
rocket_air_max_height = h(y_burn, v_burn, g, rocket_air_max_height_time, rocket_vterm)

# ================STATS================
print(f"\n=============BASIC KINEMATICS============")
print(f"Max height : {round(kinematics_max_height, 2)}m")
print(f"Flight time : {round(flight_time, 2)}s")
print(f"\n===========WITH AIR RESISTANCE===========")
print(f"Max height : {round(air_resistance_max_height, 2)}m")
print(f"Flight time : {np.round(root_air_resis, 2)}s")
print(f"Terminal velocity : {round(v_term, 2)} m/s")
print(f"\n==============ROCKET MOTION==============")
print(f"Height at engine shutoff : {np.round(y_burn, 2)}m")
print(f"Velocity at engine shutoff : {round(v_burn, 2)}m/s")
print(f"\nEngine burn time : {round(t_burn, 2)}s")
print(f"Time at liftoff : {round(t_liftoff, 2)}s")
print(f"\nMax height (w/o air resistance) : {round(rocket_kin_max_height, 2)}m")
print(f"Max height (w/ air resistance) : {round(rocket_air_max_height, 2)}m")

# Define arrays of time values to plot
kin_time = np.linspace(0, flight_time, 100) if flight_time else np.linspace(0, 5, 100)
air_time = np.linspace(0, root_air_resis, 100)
rocket_time = np.linspace(0, t_burn, 100)
rocket_free_fall_time = np.linspace(0, rocket_kin_flight_time, 100)
rocket_free_fall_real_time = rocket_free_fall_time + t_burn

# Plot equations
plt.plot(kin_time, f(y_initial, v_initial, g, kin_time), color='red', label='Kinematic Trajectory')
plt.plot(air_time, h(y_initial, v_initial, g, air_time, v_term), color='green', label='Kinematics w/ Air Resistance')
plt.plot(rocket_time, y(rocket_mass, g, exhaust_speed, exhaust_rate, rocket_time, t_liftoff), color='blue', label='Rocket Motion w/ Engine')
plt.plot(rocket_free_fall_real_time, f(y_burn, v_burn, g, rocket_free_fall_time), color='navy', linestyle='dotted', label='Engine Cutoff w/o Air Resistance')
plt.plot(rocket_free_fall_real_time, h(y_burn, v_burn, g, rocket_free_fall_time, rocket_vterm), color='purple', linestyle='dotted', label='Engine Cutoff w/ Air Resistance')
plt.axhline(y=0, color='black')
plt.legend()
plt.plot()
plt.xlabel("Time")
plt.ylabel("Height")
plt.show()
