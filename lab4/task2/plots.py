import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def exact_solution(x):
    return x**2 + 1  # Новое точное решение: y(x) = x^2 + 1

# Загружаем данные из файлов, созданных C++ программой
shooting_data = pd.read_csv('shooting.csv')
fd_data = pd.read_csv('finite_difference.csv')

# Интервал для точного решения: [0, 2]
x_fine = np.linspace(0, 2, 1000)
y_exact_fine = exact_solution(x_fine)

plt.figure(figsize=(12, 8))

# График 4: Сравнение методов
x_vals = shooting_data['x']
y_shooting = shooting_data['y']
y_fd_interp = np.interp(x_vals, fd_data['x'], fd_data['y'])

plt.plot(x_vals, y_shooting, 'r--', label='Метод стрельбы', linewidth=2)
plt.plot(x_vals, y_fd_interp, 'b-.', label='Конечно-разностный', linewidth=2)
plt.plot(x_fine, y_exact_fine, 'k-', label='Точное решение', linewidth=1, alpha=0.5)
plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Сравнение методов')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()