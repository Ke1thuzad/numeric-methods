import math

import numpy as np
import matplotlib.pyplot as plt


def xy_axes(subplot):
    subplot.axhline(0, color='k')
    subplot.axvline(0, color='k')


x = np.linspace(-10, 10, 300)

y1 = np.cos(x) + 0.25 * x - 0.5

y2 = np.arccos(0.5 - 0.25 * x)
y3 = 1 / np.sqrt(-(x * x) + 4 * x + 12)

y4 = -np.cos(x)

fig, ax = plt.subplots(2, 3)


# Простые итерации
ax[0, 0].set_title('Исходная $f(x)$')
ax[0, 0].plot(x, y1, label='cos(x) + 0.25x - 0.5')
ax[0, 0].set_ylim(-4, 4)

ax[0, 1].set_title(r'$\varphi(x)=arccos(0.5 - 0.25 * x)$')
ax[0, 1].plot(x, y2, label='arccos(0.5 - 0.25 * x)')
ax[0, 1].plot(x, x, label='x')
ax[0, 1].plot(x, y3, label=r'$\frac{1}{\sqrt{-x^2+4x+12}}$')
ax[0, 1].set_xlim(0, 3)
ax[0, 1].set_ylim(0, 3)

ax[0, 2].set_title(r'$\varphi^{\prime}(x)=\frac{1}{\sqrt{-x^2+4x+12}}$')
ax[0, 2].plot(x, y3, label=r'$\frac{1}{\sqrt{-x^2+4x+12}}$')


# Метод Ньютона
ax[1, 0].set_title(r'$f^{\prime\prime}(x)=-cos(x)$')
ax[1, 0].plot(x, y4, label='$-cos(x)$')

ax[1, 1].set_title(r'$f^{\prime\prime}(x)$ и $f(x)$')
ax[1, 1].plot(x, y1, label='$cos(x) + 0.25x - 0.5$')
ax[1, 1].plot(x, y4, label='$-cos(x)$')

ax[1, 2].set_title('Первый положительный корень (приближено)')
ax[1, 2].plot(x, y1, label='$cos(x) + 0.25x - 0.5$')
ax[1, 2].plot(x, y4, label='$-cos(x)$')
ax[1, 2].plot((1.5, 1.5), (-math.cos(1.5), math.cos(1.5) + 0.25 * 1.5 - 0.5), 'go', label='$x_{0}$')
# ax[1, 2].plot((1.5, np.cos(1.5) + 0.25 * x - 0.5), 'go')
ax[1, 2].set_xlim(1.2, 1.7)
ax[1, 2].set_ylim(-0.2, 0.1)


for axes in ax:
    for col in axes:
        col.legend(loc='upper center')
        col.grid(True)
        xy_axes(col)


plt.show()
