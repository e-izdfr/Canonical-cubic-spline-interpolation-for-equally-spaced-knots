import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sym
def canonical_spline(matrix):
    def f(x):
        return 1 / (1 + (x**2))
    n = len(matrix) - 1
    h = (matrix[-1] - matrix[0]) / n
    L = []
    for i in matrix:
        L.append(f(i))
    l = [0]
    A = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        l.append(12 * ((((L[i + 1] - L[i]) / (matrix[i + 1] - matrix[i])) - ((L[i] - L[i - 1]) / (matrix[i] - matrix[i - 1]))) / 
            (matrix[i + 1] - matrix[i - 1])))
        A[i - 1] = np.append(np.append(np.zeros(i - 1), np.array([1, 4, 1])), np.zeros(n - i - 1))
    A[n - 1] = np.append(np.array([1]), np.zeros(n))
    A[n] = np.append(np.zeros(n), np.array([1]))
    l = l + [0, 0]
    del(l[0])
    l = np.array(l)
    m = np.matmul(np.linalg.inv(A), l)
    d = []
    c = []
    for i in range(0, n):
        d.append(L[i] - ((h**2) / 6) * m[i])
        c.append(((L[i + 1] - L[i]) / h) - ((h / 6) * (m[i + 1] - m[i])))
    SUB = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
    x = sym.symbols('x')
    for i in range(0, n):
        print(f's{str(i).translate(SUB)} is : ', sym.Poly([(-1 / 6 * h) * (m[i + 1] - m[i]), (-1 / 6 * h) * (3 * m[i + 1] * matrix[i] - 3 * m[i] * matrix[i + 1]), 
        (-1 / 6 * h) * (3 * m[i] * (matrix[i + 1]) ** 2 - 3 * m[i + 1] * (matrix[i]) ** 2 + (-6) * h * c[i]), (-1 / 6 * h) * (-1 * m[i] * (matrix[i + 1]) ** 3) + 
        (matrix[i] * (matrix[i + 1]) ** 3 * m[i + 1] + (-6 * h) * (d[i] - c[i] * matrix[i]))], x).as_expr())
    def S(x):
        if matrix[n - 1] <= x:
            return ((-(x - matrix[n])**3) / (6 * h)) * m[n - 1] + (((x - matrix[n - 1])**3) /(6 * h)) * m[n] + c[n - 1] * (x - matrix[n - 1]) + d[n - 1]
        else:
            for i in range(0, n):
                if matrix[i] <= x:
                    pass
                else:
                    j = i - 1
                    return ((-(x - matrix[j + 1])**3) / (6 * h)) * m[j] + (((x - matrix[j])**3) /(6 * h)) * m[j + 1] + c[j] * (x - matrix[j]) + d[j]
    def e(x):
        return S(x) - f(x)
    e = np.vectorize(e)
    x = np.arange(matrix[0], matrix[-1], 0.00001)
    S = np.vectorize(S)
    cursor = matplotlib.widgets.Cursor(plt.figure(1).subplots(), horizOn = True, vertOn=True, color='r', linewidth=1, useblit=True)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('canonical cubic spline interpolation for equally-spaced knots', c='g')
    plt.plot(x, S(x), label='spline')
    plt.plot(x, e(x), c='g', label='error')
    plt.plot(x, f(x), c='c', label='function')
    plt.legend()
    plt.scatter(matrix, L, c='r')
    for i in range(0, n):
        plt.text(matrix[i] - 0.0001, L[i] + 0.0001, f'({matrix[i]},{L[i]})', fontsize=10)
    plt.text(matrix[n] - 0.0003, L[n] - 0.0001, f'({matrix[n]},{L[n]})', fontsize=10)
    plt.show()
canonical_spline([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])