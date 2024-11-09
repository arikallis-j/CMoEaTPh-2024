struct Equation
    name::String 
    N::Int
    M::Int 
    a::Float64 
    q::Float64
    X_rng 
    T_rng
    f_cond
    u_cond
end

function calculate(eq, num, args)
    println(eq, num)
    println(args)

    N = args["t"]
    M = args["x"]
    a = args["a"]
    q = args["q"]

    L = 1.0
    h = L / M

    T = abs(L / a)
    dt = T / N

    kurant = a * dt / h


    equation = Equation("Test", N, M, a, q, (0,L), (0,T), (phi,psi), (5, 10))

    println(equation.name)
end


# Создание структуры коэффициентов полинома
struct Derivative{R}
    coeffs::Matrix{R}
    point
end

# Объявление "вызова" для типа Polynomial
function (d::Derivative)(y)
    C = d.coeffs
    a, b = d.point[1], d.point[2]
    dy = zeros(size(y))
    A, B = size(C)
    N, M = size(y)
    @show N, M
    @show A, B 
    @show a, b
    for n=1:N
        for m=1:M
            for k = 1:A
                for l = 1:B
                    i, j = k - a,  l - b
                    if ((n+i) <= 0) || ((n + i) >= N) || ((m+j) <= 0) || ((m+j) >= M)
                        @show "continue", i, j, n, m
                        continue
                    end
                    @show i, j, n, m
                    @show (n+i), (m+j)
                    dy[n,m] += C[k,l] * y[n + i , m + j]
                    @show dy[n,m]
                end
            end
        end
    end
    return dy
end

Ct =  [-1 ; 
        1 ;;]
@show Ct
at = (1, 1)
Cx = [ -1 1;]
@show Cx
ax = (1, 2)
@show Cx[ax[1],ax[2]]
R = 3


dt = Derivative(Ct, at)
dx = Derivative(Cx, ax)

y::Matrix = ones(1, R)
for m = 1:R
    x = (m-1)
    y[1:end,m] =  (x*x + 1) .+ zeros(size(y[1:end,m]))
end
@show y

dt_y = dt(y)
dx_y = dx(y)
@show dt_y
@show dx_y

# sign(x) = (x ≥ 0) ? "+" : "-"


# Base.ndims(p::Polynomial) = length(p.coeffs) - 1

# # Задание случая отсутствия аргументов функции
# (p::Polynomial)() = p(0)

# # Функции для коэффициентов суммы полиномов
# function psum(a,b)
#     z = (length(a)>length(b)) ? copy(a) : copy(b)
#     for k=1:min(length(a), length(b))
#         z[k] = a[k] + b[k]
#     end
#     return z
# end

# # Функции для коэффициентов произведения полиномов
# function pmul(a, b)
#     len = length(a)+length(b)-1
#     a1 = psum(a, zeros(Int64, len))
#     b1 = psum(b, zeros(Int64, len))
#     z = zeros(Int64, len)
#     for k=1:len
#         for i=1:k
#             z[k] += a1[i]*b1[k-i+1]
#         end
#     end
#     return z 
# end

# # Объявление сложения, умножения и возведения в степень полиномов
# import Base: *, +, ^
# +(f::Polynomial, g::Polynomial) = Polynomial(psum(f.coeffs, g.coeffs)) # Cложение
# +(α::Number, g::Polynomial) = Polynomial(psum([α], g.coeffs)) # Cложение
# *(f::Polynomial, g::Polynomial) = Polynomial(pmul(f.coeffs, g.coeffs)) # Умножение
# *(α::Number, g::Polynomial) = Polynomial(α.*g.coeffs) # Умножение
# ^(f::Function, n::Integer) = n == 1 ? f : f*f^(n-1)

# # Объявление полиномиальной функции
# p = Polynomial([1,10,100])
# x = Polynomial([0,1])
# @show p
# @show ndims(p)
# @show p(3)
# p1 = x^2
# p2 = x - 1
# @show 2*p1 - p2 + 3
# @show p2^5;

# # Экспонента через ряд Тейлора
# pow(n) = Polynomial([ i!=(n+1) ? 0 : 1 for i=1:(n+1)])
# myexp = sum(1//factorial(big(n)) * pow(n) for n in 0:100);  
# [myexp(1); exp(big(1))];