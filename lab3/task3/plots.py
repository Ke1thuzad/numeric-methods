import numpy as np
import matplotlib.pyplot as plt

x_data = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1])
y_data = np.array([-2.3026, -0.69315, -0.10536, 0.26236, 0.53063, 0.74194])

poly1_coeffs = [-1.7745, 1.3758]

poly2_coeffs = [-2.4604, 3.4061, -0.9229]

def poly1(x):
    return poly1_coeffs[0] + poly1_coeffs[1] * x

def poly2(x):
    return poly2_coeffs[0] + poly2_coeffs[1] * x + poly2_coeffs[2] * x**2

x_plot = np.linspace(0.05, 2.15, 400)
y_poly1 = poly1(x_plot)
y_poly2 = poly2(x_plot)

plt.figure(figsize=(12, 6))

plt.plot(x_plot, y_poly1, 'r-', label='Полином 1-й степени: -1.7745 + 1.3758x')
plt.plot(x_plot, y_poly2, 'g-', label='Полином 2-й степени: -2.4604 + 3.4061x - 0.9229x²')
plt.plot(x_data, y_data, 'ko', label='Исходные точки', zorder=5)

plt.grid(True, alpha=0.3)
plt.legend()
plt.title('Аппроксимация полиномами 1-й и 2-й степени')
plt.xlabel('x')
plt.ylabel('y')

plt.tight_layout()
plt.show()