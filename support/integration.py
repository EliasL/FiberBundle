from sympy import *
import numpy as np

x_2 = symbols("x_2")
x_k = symbols("x_k")
x_m = symbols("x_m")
x_l = symbols("x_l")
t_0 = symbols("t_0")
N = symbols("N")
#l = symbols("l")
#m = symbols("m")
m=2
l=2

C = symbols("C")

def P(x):
    return  (x-t_0)/(1-t_0)

def p(x):
    return 1/(1-t_0)


def phi(x_m, x_l):
    return C * P(x_m)**(m-1) * (P(x_l)-P(x_m))**(l-m-1) * (1-P(x_l))**(N-l)*p(x_m)*p(x_l)

print("Calculating...")
integral = integrate(phi(x_2, x_k), (x_k,0,1), (x_2,2*x_k/3, x_k))
print(integral.doit())

C*Integral((t_0 - x_2)*Integral((-x_k/(1 - t_0) + 1/(1 - t_0))**N/(x_2*x_k**2 - 2*x_2*x_k + x_2 - x_k**3 + 2*x_k**2 - x_k), (x_k, 0, 1)), (x_2, 2*x_k/3, x_k))