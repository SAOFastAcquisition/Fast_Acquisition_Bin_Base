import numpy as np
import matplotlib.pyplot as plt
import allantools 

# Генерируем тестовые данные
np.random.seed(42)
n = 10000
t = np.linspace(0, 100, n)

# Создаем сигнал с разными типами шумов
white_noise = np.random.normal(0, 0.1, n)
flicker_noise = np.cumsum(np.random.normal(0, 0.005, n)) * 0.2
drift = 0.0001 * t
data = white_noise + flicker_noise + drift

# Вычисляем частоту дискретизации
sample_rate = 1 / (t[1] - t[0])

# ВЫЧИСЛЕНИЕ ВАРИАЦИИ АЛЛАНА
taus, adev, adev_error, n = allantools.adev(
    data,
    rate=sample_rate,
    taus='decade'  # автоматический выбор интервалов
)

# Построение графика
plt.figure(figsize=(12, 8))
plt.loglog(taus, adev, 'o-', linewidth=2, markersize=6)
plt.xlabel('Время усреднения, τ (сек)')
plt.ylabel('Отклонение Аллана, σ(τ)')
plt.title('Анализ Аллана с помощью allantools')
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.show()

# Дополнительная информация
print("Результаты анализа:")
print(f"Количество точек анализа: {len(taus)}")
print(f"Диапазон τ: {taus[0]:.3f} - {taus[-1]:.1f} сек")
print(f"Минимальное отклонение: {np.min(adev):.6f}")