import numpy as np
import matplotlib.pyplot as plt
import allantools

# Создаем данные с известными характеристиками
np.random.seed(42)
n = 8000
t = np.linspace(0, 80, n)

# Сигнал с преобладающим фликкер-шумом
data = np.cumsum(np.random.normal(0, 0.01, n)) * 0.1

# Настраиваем параметры анализа
sample_rate = 100  # Гц

# Ручное задание интервалов τ
custom_taus = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]

# Вычисление с кастомными τ
taus, adev, adev_error, n = allantools.adev(
    data,
    rate=sample_rate,
    taus=custom_taus  # явное задание интервалов
)

# Детальный анализ
plt.figure(figsize=(12, 8))

# График отклонения Аллана с ошибками
plt.errorbar(taus, adev, yerr=adev_error,
            fmt='o-', capsize=5, linewidth=2,
            label='ADEV с ошибками')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Время усреднения, τ (сек)')
plt.ylabel('Отклонение Аллана, σ(τ)')
plt.title('Детальный анализ стабильности')
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend()

# Добавляем аннотации
for i, (tau_val, adev_val) in enumerate(zip(taus, adev)):
    if i % 2 == 0:  # аннотируем через одну точку
        plt.annotate(f'τ={tau_val:.1f}s',
                    xy=(tau_val, adev_val),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8)

plt.show()

# Вывод статистики
print("=== СТАТИСТИКА АНАЛИЗА ===")
print(f"Общее время данных: {t[-1]:.1f} сек")
print(f"Частота дискретизации: {sample_rate} Гц")
print(f"Количество точек анализа: {len(taus)}")
print("\nРезультаты по интервалам τ:")
for tau, sigma, error in zip(taus, adev, adev_error):
    print(f"τ={tau:6.1f} сек: σ={sigma:.6f} ± {error:.6f}")