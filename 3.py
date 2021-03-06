import numpy as np
import math

n = 9
A = [
    [50. + 3 * n, 10. - n, 3.],
    [10. - n, 20. + 2 * n, 10. - n],
    [3., 10. - n, 90. - n]
]


def rotation_metod(matrix, eps=0.01):
    a = np.array(matrix)
    n = len(matrix)
    a_prev = np.zeros((n,n))
    h = np.eye(n)
    h_box = []
    v = np.zeros((n,n))
    max_el = a[0][1]
    x = 0
    y = 1
    k = 0
    while max_el > eps:
        max_el = a[0][1]
        x = 0
        y = 1
        a_prev = a
        # Находим максимальный недиагональный элемент
        for i in range(n):
            for j in range(n):
                if i != j and abs(a_prev[i][j]) > max_el:
                    x = i
                    y = j
                    max_el = a_prev[i][j]
        if x > y:
            x, y = y, x
        # Определяем угол Фи
        fi = 1/2 * math.atan((2 * a_prev[x][y]) / (a_prev[x][x] - a_prev[y][y]))
        # Заполняем матрицу H
        h = np.eye(n)
        h[x][x] = round(math.cos(fi), 5)
        h[x][y] = round(- math.sin(fi), 5)
        h[y][y] = round(math.cos(fi), 5)
        h[y][x] = round(math.sin(fi), 5)
        # Находим следующее приближение A
        a = np.matmul(np.matmul(h.T, a_prev), h)
        # Добавляем матрицу H в список
        h_box.append(h)
        k += 1
        print(f'\nA{k}:')
        print(a)
        print(f'H{k}:')
        print(h)
    sob_znach = [str(a[i][i]) for i in range(n)]
    v = np.eye(n)
    # Перемножаем все матрицы H и получаем V - матрицу с собственными векторами
    for matr in h_box:  
        temp = np.matmul(v, matr)
        v = temp
    sob_vec = []
    # Выделяем вектора из матрицы и проводим нормировку
    for i in range(n):
        vec = v[:,i]/ v[i][i]
        sob_vec.append(vec)
    print('\nСобственные значения:')
    print('\n'.join(sob_znach))
    print('\nСобственные вектора:')
    for i, vec in enumerate(sob_vec):
        print(f'x[{i}] = ',vec)

rotation_metod(A, eps=0.01)