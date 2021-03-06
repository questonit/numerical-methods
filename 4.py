from math import log, exp


# Исходная функция
def f(x):
    return (2 - x) * exp(x)

# Производная исходной функции
def df(x):
    return (2 - x) * exp(x) - exp(x)

# Вторая производная исходной функции
def d2f(x):
    return (2 - x) * exp(x) - 2 * exp(x)

# Так как при попытке представить исходную функцию в виде, подходящем для итераций, получается x = 2,
# то для метода итераций используем  другую функцию
def f1(x):
    return log(x) + 2 * x**2 - 6
# Рекуррентная функция для итераций
def g1(x):
    return x - 0.1 * f1(x)

# Производная рекуррентной функции
def dg1(x):
    return -0.4 * x - 0.1 / x + 1

def iterations(f, g, dg, start, eps):
    # Массив для итераций
    x = [start]
    it = 0
    while True:
        # Ошибка если итерации не сходятся
        if dg(x[-1]) >= 1:
            raise ValueError('Iterations do not converge, try changing g function')
        # Условия выхода
        if len(x) >= 3:
            if abs(f(x[-1])) < eps and (x[-1] - x[-2])**2 / abs(2 * x[-2] - x[-1] - x[-3]) < eps:
                break
        # Добавление новой итерации
        x.append(g(x[-1]))
        it += 1
        print('Итерация', it, ': x =', g(x[-1]))
    return x[-1]


def dichotomy(f, start, end, eps):
    # Границы поиска
    a = start
    b = end
    # Проверка на наличие рещения
    if f(a) * f(b) > 0:
        raise ValueError('No solutions. Perhaps the borders are not selected correctly.')
    it = 0
    while b - a > eps * 2:
        # Подсчет середины
        c = (a + b) / 2
        it += 1
        print('Приближение', it, ': x =', c)
        # Установка новых границ
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2


def newton(f, df, d2f, start, end, eps):
    # Выбор начального приближения
    x0 = end if f(start) * d2f(start) > 0 else start
    # Первое приближение
    x1  = x0 - f(x0)/df(x0)
    # Следующее приближение пока не достигнута точность
    it = 0
    while abs(x0-x1) > eps:
        x0 = x1
        x1 = x1 - f(x1) / df(x1)
        it += 1
        print('Приближение', it, ': x =', x1)
    return x1


print("Метод Итераций: \n")
res = iterations(f1, g1, dg1, 0.1, 0.01)
print('Решение:', res)
print("Метод Дихотомии:\n")
res = dichotomy(f, 0.1, 5, 0.01)
print('Решение:', res)
print("Метод Ньютона:\n ")
res = newton(f, df, d2f, -1, 4, 0.01)
print('Решение:', res)