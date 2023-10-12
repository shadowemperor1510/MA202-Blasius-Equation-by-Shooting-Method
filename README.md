## MA202-Blasius-Equation-by-Shooting-Method

This is a code to solve the Blausius Equation using the shooting method. Adapted from <reference>. 
Written 25 April 2023, by Aditya Gupte and Indira Gadhvi.

The given differential equation is of the third order. So we convert it to a first order ODE. This is implemented in the definition
of the function 'f'. 

We then focus on the conditions given to us. The Blausius equation is a boundary value problem. We solve it by converting
it to an initial value problem by taking a guess for the value of the derivative of f, and solve the three first order differential equations
using the 4th order Runge-Kutta method to check if the value of the function u2 at the end point thus derived lies within a tolerance of the given final condition (called the target).
This guess is called the shooting parameter 's'.  In order to find the correct s, we adefine the function F = u2(at L for given s) - target and find its zeroes using the secant method. 
Note that we approximate infinity wtih a large value L. 

Also note that an important idea in this code was to condense all three equations and solve them simultaneously as an array, so that at each new iteration, the values of u1 and u2 could be accessed to compute u3.
