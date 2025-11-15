import numpy as np
import matplotlib.pyplot as plt

x_points = np.array([0.2, 0.6, 1.0, 1.4])
y_points = np.log(x_points)
x_star = 0.8

def lagrange_poly(x):
    term1 = 4.19124 * (x - 0.6) * (x - 1.0) * (x - 1.4)
    term2 = -3.99083 * (x - 0.2) * (x - 1.0) * (x - 1.4)
    term3 = 0.87623 * (x - 0.2) * (x - 0.6) * (x - 1.0)
    return term1 + term2 + term3

def newton_poly(x):
    term1 = -1.60944
    term2 = 2.74653 * (x - 0.2)
    term3 = -1.83683 * (x - 0.2) * (x - 0.6)
    term4 = 1.07665 * (x - 0.2) * (x - 0.6) * (x - 1.0)
    return term1 + term2 + term3 + term4

x_plot = np.linspace(0.15, 1.45, 400)
y_original = np.log(x_plot)
y_lagrange = lagrange_poly(x_plot)
y_newton = newton_poly(x_plot)

plt.figure(figsize=(12, 8))

plt.plot(x_plot, y_original, 'b-', label='ln(x)')
plt.plot(x_plot, y_lagrange, 'r--', label='Полином Лагранжа')
plt.plot(x_plot, y_newton, 'g--', label='Полином Ньютона')
plt.plot(x_points, y_points, 'ko', label='Узлы интерполяции')
plt.plot(x_star, lagrange_poly(x_star), 'ms', label=f'X* = {x_star}')

plt.grid(True)
plt.legend()
plt.title('Интерполяционные полиномы Лагранжа и Ньютона')
plt.xlabel('x')
plt.ylabel('y')

plt.tight_layout()
plt.show()
