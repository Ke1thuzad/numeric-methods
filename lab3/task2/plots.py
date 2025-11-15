import numpy as np
import matplotlib.pyplot as plt

x = np.array([0.1, 0.5, 0.9, 1.3, 1.7])
f = np.array([-2.3026, -0.69315, -0.10536, 0.26236, 0.53063])
x_point = 1.5

splines = [
    {
        'interval': [0.1, 0.5],
        'a': -2.302600,
        'b': 4.672914,
        'c': 0.000000,
        'd': -4.058055
    },
    {
        'interval': [0.5, 0.9],
        'a': -0.693150,
        'b': 2.725047,
        'c': -4.869666,
        'd': 4.326839
    },
    {
        'interval': [0.9, 1.3],
        'a': -0.105360,
        'b': 0.906197,
        'c': 0.322540,
        'd': -0.724456
    },
    {
        'interval': [1.3, 1.7],
        'a': 0.262360,
        'b': 0.816490,
        'c': -0.546807,
        'd': 0.455672
    }
]

def cubic_spline(x_val, spline):
    dx = x_val - spline['interval'][0]
    return spline['a'] + spline['b'] * dx + spline['c'] * dx**2 + spline['d'] * dx**3

plt.figure(figsize=(10, 6))

for i, spline in enumerate(splines):
    x_start, x_end = spline['interval']

    x_interp = np.linspace(x_start, x_end, 100)
    y_interp = cubic_spline(x_interp, spline)
    plt.plot(x_interp, y_interp, label=f'Сплайн на [{x_start}, {x_end}]')

    x_extrap_left = np.linspace(max(0, x_start - 0.75), x_start, 50)
    y_extrap_left = cubic_spline(x_extrap_left, spline)
    plt.plot(x_extrap_left, y_extrap_left, '--', color=plt.gca().lines[-1].get_color(), alpha=0.3)

    x_extrap_right = np.linspace(x_end, min(2.2, x_end + 0.75), 50)
    y_extrap_right = cubic_spline(x_extrap_right, spline)
    plt.plot(x_extrap_right, y_extrap_right, '--', color=plt.gca().lines[-1].get_color(), alpha=0.3)

plt.plot(x, f, 'ko', markersize=6, label='Узлы интерполяции')
plt.plot(x_point, cubic_spline(x_point, splines[3]), 'ro', markersize=8, label=f'f({x_point})')

plt.grid(True)
plt.legend()
plt.title('Кубические сплайны')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0, 2.2)
plt.ylim(-2.5, 1.0)

plt.tight_layout()
plt.show()