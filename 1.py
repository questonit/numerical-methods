import numpy as np
import math


a = [
    [0.54, -0.04, 0.10],
    [-0.04, 0.50, 0.12],
    [0.10, 0.12, 0.71]
]
b = [0.33, -0.05, 0.28]

# Метод Гаусса
def gauss(left, right, prec=3):
    # Создаем расширенную матрицу
    arr = np.concatenate((np.array(left), np.array([right]).T), axis=1)
    print('\nИсходная матрица:')
    print(arr)
    # Проверка совместности
    if np.linalg.matrix_rank(left) != np.linalg.matrix_rank(arr):
        return 'Решений нет!'
    # Приводим к ступенчатому виду
    for j in range(len(arr)):
        # Находим ведущий элемент
        lead = j
        for i in range(j, len(arr)):
            if (arr[i][j] > arr[lead][j] and arr[i][j] != 0):
                lead = i
        # Если все элементы строки - 0, пропускаем итерацию
        if arr[lead][j] == 0:
            continue
        # Выносим строку с ведущим элементом вверх
        arr[[j, lead]] = arr[[lead, j]]
        # Обнуляем нижестоящие элементы
        arr[j] = arr[j] / arr[j][j]
        for i in range(j + 1, len(arr)):
            arr[i] = arr[i] - arr[j] * arr[i][j]
        print('\nШаг ', j)
        print(arr)
    # Приводим матрицу к единичной
    for j in reversed(range(len(arr))):
        for i in reversed(range(j)):
            arr[i] = arr[i] - arr[j] * arr[i][j]
    print('\nМатрица в единичном виде')
    print(arr)
    # Формируем и возвращаем результат
    answer = {('x' + str(i + 1))
               : format(arr[:, -1][i], f'.{prec}f') for i in range(len(arr))}
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
    return norm_1(matrix) < 1 or norm_2(matrix) < 1 or norm_3(matrix) < 1

# Метод простой итерации
def iteration(left, right, eps=0.0001, prec=5):
    # Формируем матрицу Альфа
    alpha = [[(-left[i][j] / left[i][i]) if (i != j)
              else 0 for j in range(len(left))] for i in range(len(left[0]))]
    # Формируем вектор Бета
    beta = np.array([right[i] / left[i][i] for i in range(len(left))])
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
            x = np.dot(alpha, prev_x) + beta
            # Считаем точность
            cur_eps = cur_eps * norm_alpha
            it += 1
            print('Итерация', it, ': X =', x)
        # Формируем и возвращаем результат
        answer = {('x' + str(i + 1))
                   : format(x[i], f'.{prec}f') for i in range(len(x))}
        return answer
    # Если решение не сходится - ошибка
    else:
        return 'Решение не сходится!'

print('Метод Гаусса')
res = gauss(a, b, prec=5)
print('Решение:', res)
print('\nМетод простой итерации')
res = iteration(a, b, eps=0.01, prec=5)
print('Решение:', res)

