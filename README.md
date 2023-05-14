# Rocket Trajectory Calculator 
This command line-based python simulator plots rocket and projectile motion
using calculus-derived equations. This is an attempt to illustrate the
difference in maximum altitude and motion behaviour when considering factors
such as air resistance and non-constant forces.

## Program Overview
### Kinematic Projectiles
Grade 11 & 12 physics teaches us about projectiles and how the maximum 
height of a projectile is dependent on the initial speed of the object. We 
can use the basic kinematic equation `yf = yi + vi * t + 0.5 * a * t^2` to
plot y vs time in a graph. This is however an ideal scenario that only
applies to a system without aerodynamic drag and with initial speed being the 
factor countering gravity's pull for a short duration of time. 

### Kinematic Projectiles with Air Resistance
We can make this model more accurate by including air resistance 
`Fa = c * v^2` as another downward force which will reduce the maximum
height that a projectile can reach. Using calculus to solve for y vs time, 
we see that air resistance only becomes a significant contributing force
when the initial velocity is close to the terminal velocity of the object.

For a given initial position and velocity, the user can choose to only plot
these two projectiles without the rocket motion part of the program. 

### Rocket Motion
While projectiles start with a given initial speed and slow down as they 
gain altitude in a parabolic arc, rockets gain speed exponentially as almost
the entirety of the rocket's mass is propellant which is lost through the
engine which expels the fuel to generate thrust. During liftoff, the rocket
gains altitude in an exponential arc until engine shutoff (when the rocket
runs out of fuel).

After engine shutoff, the rocket becomes a projectile and obeys the same
laws of motion as a normal projectile which we can also graph using the 
equations given earlier. This program will plot the rocket trajectory
after engine shutoff using the kinematic equations with and without air
resistance. 

## Usage
While many rocket trajectory calculators take thrust and impulse as inputs 
for the rocket engine, I chose to use **Exhaust Speed** and **Exhaust Rate**
which allowed me to derive the rocket motion equations by hand. 

**Exhaust Speed** is the speed at which mass is ejected from the rocket (m/s). 
- For a rocket weighing ~1-5kg, **Exhaust Speed (u)** is ~100-300m/s*
- For a rocket weighing ~4,000,000kg, **Exhaust Speed (u)** is ~2,000m/s*

**Exhaust Rate** is how much mass is ejected from the rocket per unit of
time (kg/s).
- For a rocket weighing ~1-5kg, **Exhaust Rate (alpha)** is ~0.05kg/s*
- For a rocket weighing ~4,000,000kg, **Exhaust Rate (alpha)** is ~20,000kg/s*

*Exhaust Rate and Speed cannot be arbitrary values as an engine with an 
exhaust rate too large for the amount of fuel will never achieve liftoff.
Similarly, an engine with an exhaust speed too small for the weight of the
vehicle will also never achieve liftoff. The values provided illustrate an
accurate ratio between rocket mass, exhaust speed and exhaust rate. 

**Drag Coefficient** is a unitless value that is affected by the 
aerodynamic shape of the vehicle. It ranges between 0.2 and 2.
- A cubic object has a drag coefficient of 2.08
- A spherical object has a drag coefficient of 0.38
- Rockets generally have a drag coefficient of ~0.75

## Installation
This program uses several python libraries that should be installed. 
- Install **Numpy** using `pip install numpy`
- Install **MatPlotLib** using `pip install matplotlib`
- Install **SciPy** using `pip install scipy`

## About the Math
I give credit to my Grade 12 Physics teacher who taught us how to use 
calculus to derive the equations for non-constant forces and rocket 
motion. The derivation for the rocket motion equations is included in this
repository (if the file isn't there, I haven't found the time to cleanly 
write out the formulas).

I give credit to the following websites that helped me understand the 
equations for air resistance and rocket motion.
- For in-depth explanations on the calculations for trajectories and air
resistance : [omnicalculator.com](https://www.omnicalculator.com/)
- For derivation and explanation of air resistance : [Richard Fitzpatrick](https://farside.ph.utexas.edu/teaching/336k/Newton/node29.html)
- For rocket motion theory : [LumenLearning](https://courses.lumenlearning.com/suny-osuniversityphysics/chapter/9-7-rocket-propulsion/)
- For rocket equations and a trove of useful references : [RocketMine](https://courses.lumenlearning.com/suny-osuniversityphysics/chapter/9-7-rocket-propulsion/)
- For more information on rocket thrust and impulse : [Robert A. Braeunig](http://www.braeunig.us/space/propuls.htm)
