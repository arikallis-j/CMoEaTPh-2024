# Imports

import numpy as np

# Functions

def real_sol(X, T, c, u_0, u_1):
    N = T.size
    T_full = T[-1] - T[0]
    G = np.zeros((N, 4))
    Y = np.zeros((N, 4))
    for t in range(N):
        G[t,0], Y[t,0] = u_1, X[0]
        if c>0:
            G[t,0], Y[t,0] = u_1, X[0]
            G[t,1], Y[t,1] = u_1, c * (t/N * T_full)
            G[t,2], Y[t,2] = u_0, c * (t/N * T_full)
            G[t,3], Y[t,3] = u_0, X[-1]
        else:
            G[t,0], Y[t,0] = u_0, X[0]
            G[t,1], Y[t,1] = u_0, X[-1] + c * (t/N * T_full)
            G[t,2], Y[t,2] = u_1, X[-1] + c * (t/N * T_full)
            G[t,3], Y[t,3] = u_1, X[-1]
    return G, Y

def equal(F, t, x):
    return F[t, x]

def phi(x, u_0, u_1, point=0):
    return u_1 if x==point else u_0

def psi(t, u_1):
    return u_1

def phi_m(r_m):
    if r_m>1:
        return min(2, r_m)
    elif (r_m>0 and r_m<=1):
        return min(1, 2*r_m)
    else:
        return 0
    
