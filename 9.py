import sympy
import math


def set_h(f, a, b, eps):
    x = sympy.symbols('x')
    df4 = sympy.diff(f, x, 4)
    i = a
    max_df = 0.
    while (i <= b):
        pr = df4.subs(x , i)
        max_df = max(max_df, pr)
        i += 0.1
    temp = (180 * eps) / ((b - a) * max_df)
    h = pow(temp, 1 / 4)
    return h



def simpson(f, a, b, eps):
    h = set_h(f, a, b, eps)
    # h = 0.001
    x = sympy.symbols('x')
    points = []
    n = round((b - a) / h)
    if n % 2:
        n += 1
    m = n // 2
    arg = a
    # вычислим все y и занесем точки в таблицу
    while arg <= b:
        y = f.subs(x, arg)
        points.append((arg, y))
        arg += h
    m = len(points) // 2
    first_sum = 0
    second_sum = 0
    # вычислим суммы для нечетных и четных точек
    for i in range(1, m):
        first_sum += points[2 * i - 1][1]
        second_sum += points[2 * i][1]
    first_sum += points[2 * m - 1][1]
    # подставим значения в формулу
    integr = (h / 3) * (points[0][1] + points[n-1][1] + 4 * first_sum + 2 * second_sum)
    return integr


x = sympy.symbols('x')
f = sympy.tan(pow(x, 2) + 0.5) / (1 + 2 * x)
print('Метод Симпсона: ', simpson(f, 0.5, 1., 0.001))