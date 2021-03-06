import numpy as np
import math


a = [
    [ 2., -8.,  0.,  0.],
    [-2.,  6., -1.,  0.],
    [ 0., -4., 10.,  2.],
    [ 0.,  0.,  6., -1.]
]
b = [11.4, 5.7, 15.2, 9.5]

# Система, для которой метод Зейделя сходится
# a = [
#     [10.,1.,1.],
#     [2.,10.,1.],
#     [2.,2.,10.]
# ]

# b = [12., 13., 14.]

# Метод прогонки
def sweep(left, right, prec = 3):
    for i in range(len(left)):
        for j in range(len(left[0])):
            if (i < j - 1 or i > j + 1) and left[i][j] != 0:
                return 'Матрица не трёхдиагональная!'
    # Инициализация вспомогательных массивов P[i] и Q[i]
    p = np.zeros(len(left) + 1)
    q = np.zeros(len(left) + 1)
    # Заполнение P[i] и Q[i]
    for i in range(len(left)):
        # Подсчет необходимых коэффициентов
        a = 0 if i == 0 else left[i][i-1]
        b = left[i][i]
        c = 0 if i == len(left) - 1 else left[i][i+1]
        d = right[i]
        # Вычисление P[i] и Q[i]
        p[i+1] = -c / (b + a * p[i])
        q[i+1] = (d - a * q[i]) / (b + a * p[i])
    print('\nТаблица коэффициентов P:')
    print(p)
    print('\nТаблица коэффициентов P:')
    print(q)
    # Обратный ход прогонки (нахождение решений)
    x = np.full(len(left), q[-1])
    for i in reversed(range(1, len(q) - 1)):
        x[i - 1] = q[i] + p[i] * x[i]
    # Формируем и возвращаем результат
    answer =  {('x' + str(i + 1)) : format(x[i], f'.{prec}f') for i in range(len(x))}
    return answer


def norm_1(matrix):
    data = np.array(matrix)
    return max([np.sum(np.absolute(data[i])) for i in range(len(data))])


def norm_2(matrix):
    data = np.array(matrix).T
    data = np.array(data)
    return max([np.sum(np.absolute(data[i])) for i in range(len(data))])


def norm_3(matrix):
    data = np.square(np.array(matrix).flatten())
    return math.sqrt(np.sum(data))


def converges(matrix):
    print('\nНормы матрицы:')
    print(norm_1(matrix), norm_2(matrix), norm_3(matrix))
    return norm_1(matrix) < 1 or norm_2(matrix) < 1 or norm_3(matrix) < 1


def seidel(left, right, eps = 0.0001, prec=5):
    # Формируем матрицу Альфа
    alpha = [[(-left[i][j] / left[i][i]) if (i != j) else 0 for j in range(len(left))] for i in range(len(left[0]))]
    # Формируем вектор Бета
    beta = np.array([ right[i] / left[i][i] for i in range(len(left))])
    # Задаем текущую точность
    norm_alpha = min(norm_1(alpha), norm_2(alpha), norm_3(alpha))
    norm_beta = norm_1(beta)
    cur_eps = norm_alpha / (1 - norm_alpha) * norm_beta
    
    # Если решение сходится
    if converges(alpha):
        # Выбираем за начальное приближение вектор Бэта
        x = np.copy(beta)
        it = 0
        # Выходим из цикла при достижении указанной точности
        while cur_eps > eps:
            # Запоминаем предыдущее значение
            prev_x = np.copy(x)
            # Считаем следующее приблеженное значение
            for i in range(len(alpha)):
                x[i] = np.dot(alpha[i], prev_x) + beta[i]
            # Считаем точность
            cur_eps = cur_eps * norm_alpha 
            it += 1
            print('Итерация', it, ': X =', x)
        # Формируем и возвращаем результат
        answer =  {('x' + str(i + 1)) : format(x[i], f'.{prec}f') for i in range(len(x))}
        return answer
    # Если решение не сходится - ошибка
    else:
        return 'Решение не сходится!'
        

print('Метод прогонки:')
res = sweep(a, b, prec=5)
print('\nРешение:', res)
print('Метод Зейделя:')
res = seidel(a, b, eps = 0.01, prec = 5)
print('\nРешение:', res)