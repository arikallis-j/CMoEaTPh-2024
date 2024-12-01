from .classes import *
from .func import *

class Calculator:
    def test(self):
        print(f"Test")

    def hw_1(self, N=1):

        # Constants

        U0, U1 = 5.0, 10.0
        c = 2.0
        L = 1.0
        T = abs(L / c)

        # Main
        printer = Printer(name="first", xlim=[0,L], ylim=[U0-1,U1+1])

        """
        Good solution with 0.0<sigma<1.0 (N!=M)
        """
        N, M = 300, 200
        q = 0.15 # best result
        c = 2.0
        solver = Solver(name="test", N = N, M = M, rng_x = (0,L), rng_t = (0, T))

        # GOOD
        # Scheme_1
        def scheme_1(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            return F[t,x] - sigma * (F[t,x] - F[t,x-1])

        funcs = scheme_1, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,0), False

        # GOOD
        # Scheme_2
        def scheme_2(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            return (F[t,x] + sigma * F[t+1,x-1]) / (1 + sigma)

        funcs = scheme_2, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,0), False

        # GOOD
        # Scheme_4
        def scheme_4(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            return 0.5 * ( (F[t,x-1] + F[t, x+1]) - sigma * (F[t, x+1] - F[t,x-1]) )

        funcs = scheme_4, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), False

        # GOOD
        # Scheme_5
        def scheme_5(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            y_pos = 0.5 * sigma * (1 - sigma) * (F[t,x+1] - F[t,x])
            y_neg = 0.5 * sigma * (1 - sigma) * (F[t,x] - F[t,x-1])

            return F[t,x] - sigma * (F[t,x] - F[t,x-1]) - (y_pos - y_neg)

        funcs = scheme_5, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), False

        # GOOD
        # Scheme_5_sm
        def scheme_5_sm(F, t, x, dt, dx, a, q):
            sigma = a*dt/dx
            Dpp = F[t,x+2] - F[t,x+1]
            Dp = F[t,x+1] - F[t,x]
            Dm = F[t,x] - F[t,x-1]
            Dmm = F[t,x-1] - F[t,x-2]

            Q_pos, Q_neg = 0, 0

            if (Dpp*Dp < 0) or (Dp * Dm < 0):
                Q_pos = Dp

            if (Dmm*Dm < 0) or (Dp * Dm < 0):
                Q_neg = Dm

            return scheme_5(F, t, x, dt, dx, a) + q * (Q_pos - Q_neg)

        funcs = scheme_5_sm, phi, psi, equal
        pars = ((c,q), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-2), False


        # Printing
        # U, X, T = solver.equation(funcs=funcs, pars=pars, dt=borders_t, dr=borders_r ,invert=invert)
        # G_real, Y_real = real_sol(X, T, c, U0, U1)
        # printer.animate(X, U, Y=Y_real, G=G_real)

        """
        Good solution with -1.0<sigma<0.0 (N!=M, c<0)
        """
        N, M = 300, 200
        c = -2.0
        solver = Solver(name="test", N = N, M = M, rng_x = (0,L), rng_t = (0, T))

        # GOOD
        # Scheme_3
        def scheme_3(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            return (F[t,x] - sigma * F[t+1,x+1]) / (1 - sigma)

        funcs = scheme_3, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), True

        # Printing
        # U, X, T = solver.equation(funcs=funcs, pars=pars, dt=borders_t, dr=borders_r ,invert=invert)
        # G_real, Y_real = real_sol(X, T, c, U0, U1)
        # printer.animate(X, U, Y=Y_real, G=G_real)

        """
        Good solution with sigma=1.0 (N=M)
        """
        N, M = 500, 500
        c = 2.0
        eps = 1e-15 # best result
        solver = Solver(name="test", N = N, M = M, rng_x = (0,L), rng_t = (0, T))

        # GOOD
        # Scheme_6
        def scheme_6(F, t, x, dt, dx, a):
            sigma = a*dt/dx
            return (F[t,x] - sigma * (1.5 * F[t,x] - 2 * F[t,x-1] + 0.5 * F[t,x-2]) 
                            + 0.5 * sigma**2 * (F[t,x] - 2 * F[t,x-1] + F[t,x-2]))

        funcs = scheme_6, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), False

        # GOOD
        # Scheme_7
        def scheme_7(F, t, x, dt, dx, a, eps):
            sigma = a*dt/dx
            r_m = (F[t,x] - F[t,x-1] + eps)/(F[t,x+1] - F[t,x] + eps)
            y_pos = 0.5 * phi_m(r_m) * sigma * (1 - sigma) * (F[t,x+1] - F[t,x])
            y_neg = 0.5 * phi_m(r_m) * sigma * (1 - sigma) * (F[t,x] - F[t,x-1])

            return F[t,x] - sigma * (F[t,x] - F[t,x-1]) - (y_pos - y_neg)

        funcs = scheme_7, phi, psi, equal
        pars = ((c,eps), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), False

        # GOOD
        # Scheme_8
        def scheme_8(F, t, x, dt, dx, a):
            sigma = a*dt/dx

            def S(sigma, a, b):
                return 0.5*(1+sigma)*a + 0.5*(1-sigma)*b
            def N(sigma, a, b):
                return 0.5*(3-sigma)*a - 0.5*(1-sigma)*b

            if np.abs(F[t,x+1] - F[t,x]) >= np.abs(F[t,x] - F[t,x-1]):
                y_pos = N(sigma, F[t,x], F[t,x-1])
            else:
                y_pos = S(sigma, F[t,x], F[t,x+1])

            if np.abs(F[t,x] - F[t,x-1]) >= np.abs(F[t,x-1] - F[t,x-2]):
                y_neg = N(sigma, F[t,x-1], F[t,x-2])
            else:
                y_neg = S(sigma, F[t,x-1], F[t,x])

            return F[t,x] - sigma * (y_pos - y_neg)

        funcs = scheme_8, phi, psi, equal
        pars = ((c,), (U0, U1,), (U1,))
        borders_t, borders_r, invert = (0,0), (0,-1), False

        # Printing
        # U, X, T = solver.equation(funcs=funcs, pars=pars, dt=borders_t, dr=borders_r ,invert=invert)
        # G_real, Y_real = real_sol(X, T, c, U0, U1)
        # printer.animate(X, U, Y=Y_real, G=G_real)


        """
        Unlinear
        """
        N, M = 500, 200
        c = - 0.5
        solver = Solver(name="test", N = N, M = M, rng_x = (0,L), rng_t = (0, T))

        def scheme_uu(F, t, x, dt, dx, a):
            eta = a*dt/dx
            return + eta * F[t,x]*F[t,x] + F[t,x] * (1 - eta * F[t,x-1])

        funcs = scheme_uu, phi, psi, equal
        pars = ((c,), (U0, U1, ), (U1,))
        borders_t, borders_r, invert = (0,0), (0,0), False

        # Printing
        U, X, T = solver.equation(funcs=funcs, pars=pars, dt=borders_t, dr=borders_r ,invert=invert)
        printer.animate(X, U)

