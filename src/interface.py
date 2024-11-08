from .classes import *
from .func import *

class Calculator:
    def test(self):
        print(f"Test")

    def hw_1_1(self, x=100, t=100, n=1, a=2, q=0.25):
        N = 300
        M = 200

        L = 1.0
        h = L / M

        T = abs(L / a)
        dt = T / N

        kurant = a * dt / h

        equation = Equation(name="Test", N=N, M=M, a=a, 
                            X_rng=(0,L), T_rng=(0,T), 
                            f_cond = (phi,psi), u_cond = (5, 10),
                            q = q)

        equation.conditional()
        if n==1:
            print("scheme 1")
            equation.scheme_1()
        elif n==2:
            print("scheme 2")
            equation.scheme_2()   
        elif n==3:
            print("scheme 3")
            equation.scheme_3()
        elif n==4:
            print("scheme 4")   
            equation.scheme_4()
        elif n==5:
            print("scheme 5")   
            equation.scheme_5()
        elif n==5.1:
            print("scheme 5.1")   
            equation.scheme_5_sm()
        elif n==6:
            print("scheme 6")   
            equation.scheme_6()
        elif n==7:
            print("scheme 7")   
            equation.scheme_7()
        elif n==8:
            print("scheme 8")   
            equation.scheme_8()
        equation.draw()
        # for k in range(equation.N):
        #     print(equation.Y[k,:])
    
    def hw_1_2(self, x=300, t=300, n=1, a=2, q=0.25):
        N = x
        M = t

        L = 1.0
        h = L / M

        T = abs(L / a)
        dt = T / N

        kurant = a * dt / h

        equation = Equation(name="Test", N=N, M=M, a=a, 
                            X_rng=(0,L), T_rng=(0,T), 
                            f_cond = (phi,psi), u_cond = (5, 10),
                            q = q)
        equation.conditional()
        if n==1:
            print("scheme 2-1")
            equation.scheme_2_1()

        equation.draw()
        # for k in range(equation.N):
        #     print(equation.Y[k,:])