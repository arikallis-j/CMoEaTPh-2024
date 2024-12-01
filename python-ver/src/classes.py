# Imports

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Classes

class Func:
    def __init__(self, func, pars):
        self.func = func
        self.pars = pars
    def __call__(self, *args):
        return self.func(*args, *self.pars)
    
class Setting:
    def __init__(self, standart_border=True, borders_t = (0,0), borders_r = (0,0)):
        self.standart_border = standart_border
        self.t_left = borders_t[0]
        self.t_right = borders_t[1]
        self.r_left = borders_r[0::2]
        self.r_right = borders_r[1::2]

class U_func:
    def __init__(self, main, inits, null, setting):
        self.null = null
        self.main = main
        self.init_t = inits[0]
        self.init_r = inits[1]
        self.setting = setting
    def __call__(self, F, *args):
        standart_border = self.setting.standart_border
        t_left = self.setting.t_left - 1
        t_right = self.setting.t_right
        r_left = self.setting.r_left
        r_right = self.setting.r_right
        

        t, r = args[0], args[1:len(args)]
        if t==t_left:
            return self.init_t(*r)
        
        i = 0
        for x_i in r:
            if standart_border:
                if x_i==r_left[i]:
                    return self.init_r(t)
            else:
                if x_i==r_right[i]:
                    return self.init_r(t)
            i+=1

        if (t>=t_right):
            return self.null(F, t, *r)
        
        i = 0
        for x_i in r:
            if (x_i>=r_right[i]):
                return self.null(F, t, *r)
            i+=1

        for x_i in r:
            return self.main(F, t, *r)

class Solver:
    def __init__(self, N = 500, M = 50, rng_x = (0,1.0), rng_t = (0, 1.0), name="name"):
        self.name = name
        self.N, self.M = N, M
        self.rng_x = rng_x
        self.rng_t = rng_t

        # spatial domain
        xmin, xmax = rng_x
        tmin, tmax = rng_t
        self.xmin, self.xmax = xmin, xmax
        self.tmin, self.tmax = tmin, tmax

        # x grid of n points
        self.X, self.dx = np.linspace(xmin,xmax,M,retstep=True)
        self.T, self.dt = np.linspace(tmin,tmax,N,retstep=True)

        # each value of the U array contains the solution for all x values at each timestep
        self.U = np.zeros((N, M))

    # explicit euler solution
    def equation(self, funcs, pars, dt=(0,0), dr=(0,0), invert=False):
        u_main, init_t, init_r, equal = funcs
        par_main, par_init_t, par_init_r = pars

        init_t_func = Func(init_t, par_init_t)
        init_r_func = Func(init_r, par_init_r)
        init_func = (init_t_func, init_r_func)
        equal_func = Func(equal, ())

        par_main = self.dt, self.dx, *par_main
        main_func = Func(u_main, par_main)

        dt_left, dt_right = dt[0], dt[1]
        dr_left, dr_right = dr[0], dr[1]

        setting = Setting(standart_border=not(invert), 
                          borders_t=(0+dt_left,self.N + dt_right), 
                          borders_r=(0+dr_left,self.M + dr_right))

        u_scheme = U_func(main_func, init_func, equal_func, setting=setting)
        for t in range(self.N):
            if invert:
                rng_r = range(self.M-1,-1,-1)
            else:
                rng_r = range(self.M)

            for r in rng_r:
                self.U[t,r] = u_scheme(self.U, t-1, r)
        return self.U, self.X, self.T

class Printer:
    def __init__(self, name="name", xlim=[0,1], ylim=[-2,2]):
        self.name = name
        self.xlim = xlim
        self.ylim = ylim
        
    def animate(self, X, U, Y=np.array([]), G=np.array([])):
        # plot solution
        plt.style.use('dark_background')
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)

        # animate the time data
        def animate_graph(i):
            ax1.clear()
            if G.size!=0:
                x0 = G[i]
                X0 = Y[i]
                plt.plot(X0,x0,color='deepskyblue') #or 'blue' or 'aqua'

            x = U[i]
            plt.plot(X,x,color='lime')

            plt.grid(True)
            plt.ylim(self.ylim)
            plt.xlim(self.xlim)

        anim = animation.FuncAnimation(fig,animate_graph,frames=len(U),interval=20)
        plt.show()
