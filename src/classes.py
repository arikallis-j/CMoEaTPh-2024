import matplotlib.pyplot as plt
import matplotlib.animation as ani

from rich.progress import track
from matplotlib import rc
from .func import *


# class System:
#     def __init__(self, equation, grid, pars):
#         self.equation = equation
#         self.grid = grid
#         self.pars = pars
#     def calc(self):
#         pass

# class Equation:
#     def __init__(self, function, conditional, rep = "F(x,t) = 0"):
#         self.F = function
#         self.cond = conditional
#         self.rep = rep
#     def __str__(self):
#         return self.rep

# class Grid:
#     def __init__(self, num, dim, field, scheme):
#         self.dim = dim
#         self.field = field
#         self.V = np.zeros(dim)
#         self.F = np.zeros(dim)
#         self.DF = np.zeros((num,len(dim),*dim))

class Equation:
    def __init__(self, name, N=3, M=5, a = 2, 
                 X_rng=(0,1), T_rng=(0,0.5), f_cond = (phi,psi), u_cond = (5, 10),
                 q = 0.25):
        self.name = name
        self.N, self.M, self.a = N, M, a
        self.L = X_rng[1] - X_rng[0]
        self.Z = T_rng[1] - T_rng[0]
        self.h, self.tau = self.L/self.M, self.Z/self.N
        self.eta = self.tau/self.h
        self.sigma = a * self.tau/self.h
        self.phi, self.psi = f_cond
        self.u_cond = u_cond
        self.q = q

        self.X_rng, self.T_rng = X_rng, T_rng
        self.I = np.full((M, N), np.arange(0, N)).T
        self.T = np.full((M, N), np.linspace(*T_rng, N)).T

        self.J = np.full((N, M), np.arange(0, M))
        self.X = np.full((N, M), np.linspace(*X_rng, M))

        self.Y = np.full((N, M), 0.0)
        self.Y0 = np.full((N, M), 0.0)

    def conditional(self):
        I, J, M = self.I, self.J, self.M
        if self.a>0:
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
            self.Y[I, J] = np.where(J==0, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
            self.x_rng = 1, M, 1
        elif self.a<0:
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
            self.Y[I, J] = np.where(J==M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
            self.x_rng = (M-1, 0, -1)

    def scheme_1(self):
        self.name = 'Явный левый уголок'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_1(i-1, j, Y, sigma), Y[i,j])
                

    def scheme_2(self):
        self.name = 'Неявный левый уголок'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_2(i-1, j, Y, sigma), Y[i,j])
    
    def scheme_3(self):
        self.name = 'Неявный правый уголок'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_3(i-1, j, Y, sigma), Y[i,j])
        
    def scheme_4(self):
        self.name = 'Схема Лакса'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            b -= 1
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
        else:
            a -= 1
            self.Y[I, J] = np.where(J==self.M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(True, scheme_4(i-1, j, Y, sigma), Y[i,j])
                
    def scheme_5(self):
        self.name = 'Схема Лакса-Вендрофа'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            b -= 1
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
        else:
            a -= 1
            self.Y[I, J] = np.where(J==self.M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_5(i-1, j, Y, sigma), Y[i,j])        
    
    def scheme_5_sm(self):
        self.name = 'Схема Лакса-Вендрофа, сглаженная'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            b -= 1
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
        else:
            a -= 1
            self.Y[I, J] = np.where(J==self.M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_5(i-1, j, Y, sigma), Y[i,j])  
        
        Y, q = self.Y, self.q
        if s==1:
            b -= 1
        else:
            a -= 1
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_5_sm(i-1, j, Y, q=q), Y[i, j])  
        
    def scheme_6(self):
        self.name = 'Схема Бима-Уорминга'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            a += 1
        else:
            b += 1
        for i in range(1,self.N-1):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(True, scheme_5(i-1, j-1, Y, sigma), Y[i,j])        
    

    def scheme_7(self):
        self.name = 'Схема TVD'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            b -= 1
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
        else:
            a -= 1
            self.Y[I, J] = np.where(J==self.M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_7(i-1, j, Y, sigma), Y[i,j])    

    def scheme_8(self):
        self.name = 'Схема ENO'
        I, J = self.I, self.J
        Y, sigma = self.Y, self.sigma
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        if s==1:
            b -= 1
            self.Y[I, J] = np.where(I==0, self.phi(self.X[I,J], u_cond=self.u_cond, point=self.T_rng[1]), self.Y[I,J])
        else:
            a -= 1
            self.Y[I, J] = np.where(J==self.M-1, self.psi(self.T[I,J], u_cond=self.u_cond), self.Y[I,J])
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_8(i-1, j, Y, sigma), Y[i,j])


    def scheme_2_1(self):
        self.name = 'Явный левый уголок, нелин.'
        I, J = self.I, self.J
        Y, eta = self.Y, self.eta
        cond = np.logical_and(I>0, J>0)
        a, b, s = self.x_rng
        for i in range(1,self.N):
            for j in range(a, b, s):
                self.Y[i, j] = np.where(i>0 and j>0, scheme_2_1(i-1, j, Y, eta), Y[i,j])

    def draw(self):
        u0, u1 = self.u_cond
        title = self.name
        fig = plt.figure()
        ax = plt.subplot()
        ln, = ax.plot([], [],)
        ln0, = ax.plot([], [],)
        ax.set_xlim(*self.X_rng)
        ax.set_ylim((0, 10 + 2))
        download = range(self.N)#track(range(self.N), description="Time: ", get_time=False)

        def update(i):
            ln.set_data(self.X[i,:], self.Y[i,:])

            x0 = [0, abs(self.a) * i / self.N * self.Z , abs(self.a) * i / self.N * self.Z , self.L]
            if self.a<0:
                x0 = self.L - np.array(x0)
            y = [u1, u1, u0, u0]
            ln0.set_data(x0, y)
            t = round(self.Z, 3)
            ax.title.set_text(title + f'\n T = {t}')
            return ln, ln0

        anime = ani.FuncAnimation(fig, update, interval=20, frames=download, blit=True)
        plt.show()
