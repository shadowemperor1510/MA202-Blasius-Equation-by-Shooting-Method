''' This is a code to solve the Blausius Equation using the shooting method. Adapted from <reference>. 
Written 25 April 2023, by Aditya Gupte and Indira Gadhvi.

The given differential equation is of the third order. So we convert it to a first order ODE. This is implemented in the definition
of the function 'f'. 

We then focus on the conditions given to us. The Blausius equation is a boundary value problem. We solve it by converting
it to an initial value problem by taking a guess for the value of the derivative of f, and solve the three first order differential equations
using the 4th order Runge-Kutta method to check if the value of the function u2 at the end point thus derived lies within a tolerance of the given final condition (called the target).
This guess is called the shooting parameter 's'.  In order to find the correct s, we adefine the function F = u2(at L for given s) - target and find its zeroes using the secant method. 
Note that we approximate infinity wtih a large value L. 

Also note that an important idea in this code was to condense all three equations and solve them simultaneously as an array, so that at each new iteration, the values of u1 and u2 could be accessed to compute u3. '''



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# We use Runge kutta order four for solving the system
def RK4_Method(f, a, L, u, N):

    """
    f : RHS fucntion of the system of ODE's
    u : the initial condition
    N : no. of interval points

    eta ranges between the domain interval: [a, L]
    
    """
    h = (L - a) / N         # step size between interval points
    answer = []             # To store the numerical solution
    interval_pts  = []      # To store the interval points
    answer.append(u)
    interval_pts.append(a)
    x = a
    for i in range(1, N + 1):       
        k1 = f(u , x)
        k2 = f(u + 0.5*h*k1  , x+0.5*h)
        k3 = f(u + 0.5*h*k2  , x+0.5*h)
        k4 = f(u + h*k3 , x+h)

        u = u + h * (k1 + 2 * k2 + 2 * k3 + k4)/6   # Value at the next point
        x = x + h                                   # The next point after step increment
        answer.append(u)
        interval_pts.append(x)
    return np.array(interval_pts ), np.array(answer)

# initial conditions with the shooting parameter s are assigned through an array 
def U_fixed(s):
    return np.array([0, 0, s])

# F(s) = u_2(s, L) - 1: Aim for s such that F(s) becomes zero.
def F(s, L):
    u = U_fixed(s)
    interval_pts , answer = RK4_Method(f, a, L, u, N)
    F = (answer[len(answer)-2])[1] - target  #this line extracts the final value at the end of the RK4 method.
    return F

#formulate as a first-order system and  define the right side function for the nonlinear ode
def f(u , t):
    f = np.zeros(3)
    f[0] = u[1]
    f[1] = u[2]
    f[2] = -0.5 * u[0] * u[2]
    return f

# the secant method for solving the algebric equation F(s) = 0
# We apply the secant method to find the zeroes for F(s,L)

''' The secant method is used to solve for the roots because it is harder to compute the value 
of the derivative of F with respect to s.'''

def secant_method(F,s0,s1,etol,L,maxiter):
  global s
  '''s0 and s1 are initial guesses.'''
  iterations = 0
  F0 = F(s0,L)  
  F1 = F(s1, L)
  while abs(F1) > etol and iterations < maxiter:
    slope = (F1 - F0)/(s1 - s0)
    s0 = s1
    s1 = s1 - F1/slope
    iterations += 1
    F0 = F1 
    F1 = F(s1,L)
  if abs(F1) > etol :
    print('max iterations reached.')
  return s1, iterations

# Parameters
a = 0
L = 10
N = 1000
target = 1

# initializing the parameter s, and setting the tolerance error
s0 = 1.0
s1 = 2.0
tol = 1.0e-6
s, tot_iters = secant_method(F, s0, s1, tol, L,200)
print('shooting parameter =', s)
print('number of iterations = ', tot_iters)
print()

# Now we use the approperaite shooting parameter in the initial condition and solve for the system of ODE
u = U_fixed(s)
interval_pts , answer = RK4_Method(f, a, L, u, N)

print("Data for f: ")
df_f = pd.DataFrame(answer[:, 0])
print(df_f,"\n")

print("Data for f': ")
df_f1 = pd.DataFrame(answer[:, 1])
print(df_f1,"\n")

print("Data for f'': ")
df_f2= pd.DataFrame(answer[:, 2])
print(df_f2,"\n")

# Ploting the results

plt.plot(interval_pts , answer[:, 0],  "r") # Plotting f
plt.plot(interval_pts ,  answer[:, 1], "g") # Plotting f'
plt.plot(interval_pts , answer[:, 2], "b")  # Plotting f"
plt.title(" Blasius Boundary layer Equation Plot  for f, f', f''")
plt.xlabel('$\eta$')
plt.ylabel("f")
plt.ylim([-3,5])
plt.legend(['f', "f'", "f''"]) # To identify the curves

plt.show()