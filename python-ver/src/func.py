import numpy as np

def scheme_1(n, m, y, sigma):
    return y[n,m] - sigma * (y[n,m] - y[n,m-1])

def scheme_2(n, m, y, sigma):
    return (y[n,m] + sigma * y[n+1,m-1]) / (1 + sigma)

def scheme_3(n, m, y, sigma):
    return (y[n,m] - sigma * y[n+1,m+1]) / (1 - sigma)

def scheme_4(n, m, y, sigma):
    return 0.5 * ( (y[n,m-1] + y[n, m+1]) - sigma * (y[n, m+1] - y[n,m-1]) )

def scheme_5(n, m, y, sigma):
    y_pos = 0.5 * sigma * (1 - sigma) * (y[n,m+1] - y[n,m])
    y_neg = 0.5 * sigma * (1 - sigma) * (y[n,m] - y[n,m-1])

    return y[n,m] - sigma * (y[n,m] - y[n,m-1]) - (y_pos - y_neg)

def scheme_5_sm(n, m, y, q):
    Dpp = y[n+1,m+2] - y[n+1,m+1]
    Dp = y[n+1,m+1] - y[n+1,m]
    Dm = y[n+1,m] - y[n+1,m-1]
    Dmm = y[n+1,m-1] - y[n+1,m-2]

    Q_pos, Q_neg = 0, 0

    if (Dpp*Dp < 0) or (Dp * Dm < 0):
        Q_pos = Dp

    if (Dmm*Dm < 0) or (Dp * Dm < 0):
        Q_neg = Dm

    return y[n+1,m] +  q * (Q_pos - Q_neg)

def scheme_6(n, m, y, sigma):
    return (y[n,m+1] - sigma * (1.5 * y[n,m+1] - 2 * y[n,m] + 0.5 * y[n,m-1]) 
                     + 0.5 * sigma**2 * (y[n,m+1] - 2 * y[n,m] + y[n,m-1]))

def scheme_7(n, m, y, sigma, eps=1e-15):
    r_m = (y[n,m] - y[n,m-1] + eps)/(y[n,m+1] - y[n,m] + eps)

    y_pos = 0.5 * phi_m(r_m) * sigma * (1 - sigma) * (y[n,m+1] - y[n,m])
    y_neg = 0.5 * phi_m(r_m) * sigma * (1 - sigma) * (y[n,m] - y[n,m-1])

    return y[n,m] - sigma * (y[n,m] - y[n,m-1]) - (y_pos - y_neg)

def scheme_8(n, m, y, sigma):
    def S(sigma, a, b):
        return 0.5*(1+sigma)*a + 0.5*(1-sigma)*b
    def N(sigma, a, b):
        return 0.5*(3-sigma)*a - 0.5*(1-sigma)*b

    if np.abs(y[n,m+1] - y[n,m]) >= np.abs(y[n,m] - y[n,m-1]):
        y_pos = N(sigma, y[n,m], y[n,m-1])
    else:
        y_pos = S(sigma, y[n,m], y[n,m+1])

    if np.abs(y[n,m] - y[n,m-1]) >= np.abs(y[n,m-1] - y[n,m-2]):
        y_neg = N(sigma, y[n,m-1], y[n,m-2])
    else:
        y_neg = S(sigma, y[n,m-1], y[n,m])

    return y[n,m] - sigma * (y_pos - y_neg)

def scheme_2_1(n, m, y, eta):
    return eta * (y[n,m])**2 + y[n,m] * (1 - eta * y[n,m-1])


def phi_m(r_m):
    if r_m>1:
        return min(2, r_m)
    elif (r_m>0 and r_m<=1):
        return min(1, 2*r_m)
    else:
        return 0

def phi(X, u_cond, point=0.0):
    u_0, u_1 = u_cond
    return np.where(X==point, u_1, u_0)

def psi(T, u_cond):
    u_0, u_1 = u_cond
    U1 = np.full(T.shape, u_1)
    return U1



 #def test_eq(I, F, DF, A):
#     """
#     equation: 
#         dF/dt + a * dF/dx = 0
#     """
#     dF = DF[0,:,:,:]
#     dFx, dFt = dF
#     a = A
#     return dFt + a * dFx

# def test_cond(V0, F, G0, A):
#     Vx0, Vt0 = V0
#     Gx0, Gt0 = G0

#     F[Vt0] =  Gt0(Vt0, A)
#     F[Vx0] =  Gx0(Vx0, A)

#     return F

# def test_scheme(V,i,j):
#     Vx, Vt = V
#     return 0

