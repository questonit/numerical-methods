import math

def f(x, y):
    return 1 + 0.6 * math.sin(x) - 1.25 * pow(y, 2)


def runge_kutta(f, a, b, x0, y0, h):
    x = x0
    y = y0
    # создадим таблицу значений
    points = []
    # добавим нулевые значения
    points.append((x, y))
    while x + h <= b:
        # вычисляем k1, k2, k3, k4 по формулам
        k1 = h * f(x, y)
        k2 = h * f(x + 1 / 2 * h, y + 1 / 2 * k1)
        k3 = h * f(x + 1 / 2 * h, y + 1 / 2 * k2)
        k4 = h * f(x + h, y + k3)
        # на их основе вычисляем дельта y
        dy = 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        # вычисляем текущий y
        x += h
        y = y + dy
        # добавляем в таблицу текущие значения
        points.append((x, y))
    return points


def adams(f, a, b, x0, y0, h, eps):
    x = x0
    y = y0
    a1 = x
    b1 = x + 3 * h
    # вычислим три узла, используя метод Рунге-Кутты
    points = runge_kutta(f, a1, b1, x, y, h)
    x, y = points[3]
    n = 3
    while x + h <= b:
        # вычислим f от четырех предыдущих точек
        f1 = f(points[n][0], points[n][1])
        f2 = f(points[n - 1][0], points[n - 1][1])
        f3 = f(points[n - 2][0], points[n - 2][1])
        f4 = f(points[n - 3][0], points[n - 3][1])
        # увеличим x
        x += h
        # вычислим y по формуле
        y = y + h / 24 * (55 * f1 - 59 * f2 + 37 * f3 - 9 * f4)
        points.append((x, y))
        # увеличим номер точки
        n += 1
    return points


x0 = 0
y0 = 0
a = x0
b = 1.
h = 0.1
eps = 0.01
result = runge_kutta(f, a, b, x0, y0, h)
print("Метод Рунге-Кутты четвертого порядка: ")
print("   x   |   y")
for (x, y) in result:
    print("%.5f" % x, "%.5f" % y)
result1 = adams(f, a, b, x0, y0, h, eps)
print()
print("Метод Адамса: ")
print("   x   |   y")
for (x, y) in result1:
    print("%.5f" % x, "%.5f" % y)