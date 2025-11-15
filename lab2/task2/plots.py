import math

import numpy as np
import matplotlib.pyplot as plt


def xy_axes(subplot):
    subplot.axhline(0, color='k')
    subplot.axvline(0, color='k')


# x1 - cosx2 = 1
# x2 - lg(x1 + 1) = 2

# x2 = arccos(x1 - 1)
# x2 = lg(x1 + 1) + 2

x1 = np.linspace(-10, 10, 2000)

x2_1 = np.acos(x1 - 1)

x2_2 = np.log10(x1 + 1) + 2


x2_der1 = np.abs(np.sin(x1))

x2_der2 = np.abs(1/((x1 + 1) * np.log(10)))

# sinx2 + 1 = 0
# x2 = arcsin(-1)

fig, ax = plt.subplots(3, 2)


# Простые итерации
ax[0, 0].set_title('Исходная система')
ax[0, 0].plot(x1, x2_1, label='$x_1 - cos(x_2) = 1$')
ax[0, 0].plot(x1, x2_2, label='$x_2 - lg(x_1 + 1) = 2$')
ax[0, 0].set_xlim(-1, 2)
ax[0, 0].set_ylim(0, 4)

ax[0, 1].set_title('Исходная система (приближен квадрат)')
ax[0, 1].plot(x1, x2_1, label='$x_1 - cos(x_2) = 1$')
ax[0, 1].plot(x1, x2_2, label='$x_2 - lg(x_1 + 1) = 2$')
ax[0, 1].set_xlim(0.4, 0.5)
ax[0, 1].set_ylim(2.1, 2.2)


ax[1, 0].set_title(r'Max суммы модулей производных $\varphi_1(x)=1+cos(x_2)$')
ax[1, 0].plot(x1, x2_der1, label='$|sin(x_2)| + 0$')


ax[1, 1].set_title(r'$\varphi_1^{\prime}(x)=|sin(x_2)|$  с коэффициентом $q=0.9$')
ax[1, 1].plot(x1, x2_der1, label='$|sin(x_2)| + 0$')
ax[1, 1].axhline(0.9, color='r', label='$q=0.9$')
ax[1, 1].axvline(2.1, linestyle='dashed', color='g', label=r'$2.1\leq x_2\leq 2.2$')
ax[1, 1].axvline(2.2, linestyle='dashed', color='g')
ax[1, 1].set_xlim(2, 2.5)
ax[1, 1].set_ylim(-0.3, 1)


ax[2, 0].set_title(r'Max суммы модулей производных $\varphi_2(x)=2+lg(x_1+1)$')
ax[2, 0].plot(x1, x2_der2, label=r'$|\frac{1}{(x_1+1)ln10}| + 0$')


ax[2, 1].set_title(r'$\varphi_2^{\prime}(x)=\frac{1}{(x_1+1)ln10}$ с коэффициентом $q=0.9$')
ax[2, 1].plot(x1, x2_der2, label=r'$|\frac{1}{(x_1+1)ln10}| + 0$')
ax[2, 1].axhline(0.9, color='r', label='$q=0.9$')
ax[2, 1].axvline(0.4, linestyle='dashed', color='g', label=r'$0.4\leq x_1\leq 0.5$')
ax[2, 1].axvline(0.5, linestyle='dashed', color='g')
ax[2, 1].set_xlim(0.2, 0.7)
ax[2, 1].set_ylim(-0.3, 1)


for axes in ax:
    for col in axes:
        col.legend(loc='best')
        col.grid(True)
        xy_axes(col)


plt.show()
