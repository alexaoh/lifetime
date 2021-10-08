import matplotlib.pyplot as plt
import numpy as np

infty = 30
x = np.linspace(0, infty, 2000)

l1 = 0.15
l2 = 0.1
l3 = 0.2
p1 = 7
p2 = 13
p3 = infty

def S(x):
    """Survival function."""
    conditions = [(0<=x) & (x < p1), (p1<=x) & (x<p2), (p2<=x)]
    functions = [lambda x: np.exp(-(x*l1)), 
                 lambda x: np.exp(-(p1*(l1-l2) + x*l2)), 
                 lambda x: np.exp(-(p1*(l1-l2) + p2*(l2-l3) + x*l3))
                ]
    return np.piecewise(x, conditions, functions)

def f(x):
    """Density function."""
    conditions = [(0<=x) & (x < p1), (p1<=x) & (x<p2), (p2<=x)]
    functions = [lambda x: l1*np.exp(-(x*l1)), 
                 lambda x: l2*np.exp(-(p1*(l1-l2) + x*l2)), 
                 lambda x: l3*np.exp(-(p1*(l1-l2) + p2*(l2-l3) + x*l3))
                ]
    return np.piecewise(x, conditions, functions)

def h(x):
    """Hazard function."""
    conditions = [(0<=x) & (x < p1), (p1<=x) & (x<p2), (p2<=x)]
    functions = [l1, l2, l3]
    return np.piecewise(x, conditions, functions)

plt.plot(x, h(x), "b")
plt.title("Hazard Function")
plt.xlabel("t")
plt.ylabel(r"$\lambda(t)$")
plt.savefig("hazard.pdf")
plt.show()

plt.plot(x, S(x), "b")
plt.title("Survival Function")
plt.xlabel("t")
plt.ylabel(r"$S(t)$")
plt.savefig("survival.pdf")
plt.show()

plt.plot(x, f(x), "b")
plt.title("Density Function")
plt.xlabel("t")
plt.ylabel(r"$f(t)$")
plt.savefig("density.pdf")
plt.show()
