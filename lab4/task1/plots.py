import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Загрузка данных
euler_data = pd.read_csv('euler.csv')
euler_cauchy_data = pd.read_csv('euler_cauchy.csv')
euler_improved_data = pd.read_csv('euler_improved.csv')
rk4_data = pd.read_csv('runge_kutta.csv')
adams_data = pd.read_csv('adams.csv')
exact_data = pd.read_csv('exact.csv')

plt.figure(figsize=(14, 8))

plt.plot(exact_data['x'], exact_data['y_exact'], 'k-', label='Точное решение', linewidth=2)
plt.plot(euler_data['x'], euler_data['y'], 'r--', label='Метод Эйлера', alpha=0.7, marker='o', markersize=4)
plt.plot(euler_cauchy_data['x'], euler_cauchy_data['y'], 'g-.', label='Метод Эйлера-Коши', alpha=0.7, marker='s', markersize=4)
plt.plot(euler_improved_data['x'], euler_improved_data['y'], 'm:', label='Улучшенный Эйлер', alpha=0.7, marker='^', markersize=4)
plt.plot(rk4_data['x'], rk4_data['y'], 'b-.', label='Рунге-Кутта 4', alpha=0.7, marker='d', markersize=4)
plt.plot(adams_data['x'], adams_data['y'], 'c:', label='Адамс 4', alpha=0.7, marker='v', markersize=4)
plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Сравнение численных методов')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()