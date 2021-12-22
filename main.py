import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.special as sc
#границы
a = 0
b = 5
p = 0
q = 5
#число точек разбиения
n = 1000
m = 1000
#шаги разбиения
hx = (b - a) / n
hksi = (q - p) / m

alpha = 1

# Функция входного сигнала
def f(x):
    return np.exp(1j * x / 10)

# Функция ядра
def K(ksi, x):
    arg = alpha*ksi*x
    b = sc.jv(2,arg)#ф-ия Бесселя 1р.2п.
    return b*x

#Отрисовка графиков
def plot_amplitude_and_phase(x, ampl, phase, fun_name, fun_arg):
    plt.plot(x, ampl)
    plt.title('Амплитуда ' + fun_name + '(' + fun_arg + ')')
    plt.xlabel(fun_arg)
    plt.ylabel('Амплитуда')
    plt.grid()
    plt.show()

    plt.plot(x, phase)
    plt.title('Фаза ' + fun_name + '(' + fun_arg + ')')
    plt.xlabel(fun_arg)
    plt.ylabel('Фаза')
    plt.grid()
    plt.show()

# Вектор f - зн-ия интегр. ф-ии в т. разбиения
fvect = np.full([n, 1], complex(0, 0))
for k in range(n):
    xk = a + k * hx
    fvect[k] = f(xk)

# Матрица А
A = np.full([m, n], complex(0, 0))
for l in range(m):
    ksil = p + l * hksi
    for k in range(n):
        xk = a + k * hx
        A[l, k] = K(ksil, xk)

# Матричная форма интегрального преобразования
F = A.dot(fvect) * hx

# Расчёт амплитуд и фаз
#входной сигнал
amplf = []
phasef = []
#интегр.преобр.
amplF = []
phaseF = []

massX = []
# входной сигнал
for k in range(n):
    xk = a + k * hx
    massX.append(xk)
    ampf = abs(f(xk))
    phf = math.atan2(f(xk).imag, f(xk).real)
    #массивы амплитуд и фаз
    amplf.append(ampf)
    phasef.append(phf)

massKsi = []
#нтегр.преоб-ие
for l in range(m):
    ksil = p + l * hksi
    massKsi.append(ksil)
    ampF = abs(F[l])
    phF = math.atan2(F[l].imag, F[l].real)
    amplF.append(ampF)
    phaseF.append(phF)


plot_amplitude_and_phase(massX,  amplf, phasef, 'f', 'x')

plot_amplitude_and_phase(massKsi,  amplF, phaseF, 'интегрального преобразования F', 'ξ')
