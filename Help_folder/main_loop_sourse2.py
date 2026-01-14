import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt

# ============================================================================
# 1. ВАШИ ЭКСПЕРИМЕНТАЛЬНЫЕ ДАННЫЕ (замените на свои!)
# ============================================================================
# Формат: массив частот в МГц и соответствующий HPBW_H в угловых секундах
# Примерные данные с вашего графика:
current_primary_file = '2025-12-21_01'

hpbw_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                           f'{current_primary_file}_hpbw.bin')
freq_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                           f'{current_primary_file}_freq.bin')
with open(hpbw_file_path, 'rb') as file:
    hpbw_h_data = pickle.load(file)
with open(freq_file_path, 'rb') as file:
    freq_data = pickle.load(file)
# freq_data = np.array([1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])  # МГц
# hpbw_h_data = np.array([200, 175, 155, 140, 128, 118, 110, 103, 97, 92, 88])  # угл. секунд

print("ВАШИ ЭКСПЕРИМЕНТАЛЬНЫЕ ДАННЫЕ:")
print("Частота (МГц) | HPBW_H (\")")
print("-" * 30)
for f, h in zip(freq_data, hpbw_h_data):
    print(f"{f:12.0f} | {h:8.1f}")

# ============================================================================
# 2. КОНСТАНТЫ И ПАРАМЕТРЫ
# ============================================================================
T_b = 20000.0  # K, яркостная температура активной области
R_s_arcsec = 150.0  # радиус источника в угловых секундах
HPBW_V0 = 75.0 * 60  # 75' в секундах на опорной частоте f0 = 1000 МГц


# ============================================================================
# 3. РАСЧЁТ ТОЛЬКО В ТОЧКАХ ИЗМЕРЕНИЙ
# ============================================================================
def calculate_T_a_direct(freq, hpbw_h):
    """
    Прямой расчёт антенной температуры для заданной частоты и HPBW_H
    Без интерполяции, используем только предоставленные значения
    """
    # HPBW по вертикали: ∝ 1/f
    hpbw_v = HPBW_V0 * (1000.0 / freq)

    # Преобразование HPBW в sigma для гауссова луча
    def hpbm_to_sigma(hpbw):
        return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    sigma_h = hpbm_to_sigma(hpbw_h)
    sigma_v = hpbm_to_sigma(hpbw_v)

    # Отношение размера источника к ширине луча
    x = R_s_arcsec / sigma_h

    # Вычисление интеграла I_src (аппроксимация)
    if x < 0.5:
        # Источник мал по сравнению с лучом
        I_src = np.pi * R_s_arcsec ** 2
    elif x > 3:
        # Луч мал по сравнению с источником
        I_src = 2.0 * R_s_arcsec * sigma_h * np.sqrt(2.0 * np.pi)
    else:
        # Промежуточный случай
        # Аппроксимация интеграла через функцию ошибок
        I_src = (np.sqrt(2 * np.pi) * sigma_h * R_s_arcsec *
                 (1.0 - np.exp(-x ** 2 / 2)) * (1.0 + 0.5 / x ** 2))

    # Площадь луча антенны
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # Антенная температура
    T_a = T_b * I_src / Omega_A

    return {
        'freq': freq,
        'hpbw_h': hpbw_h,
        'hpbw_v': hpbw_v,
        'sigma_h': sigma_h,
        'sigma_v': sigma_v,
        'R_s_sigma_h_ratio': x,
        'I_src': I_src,
        'Omega_A': Omega_A,
        'T_a': T_a
    }


# Расчёт для всех точек измерений
results = []
for i in range(len(freq_data)):
    result = calculate_T_a_direct(freq_data[i], hpbw_h_data[i])
    results.append(result)

# ============================================================================
# 4. ВЫВОД РЕЗУЛЬТАТОВ В ТАБЛИЦЕ
# ============================================================================
print("\n" + "=" * 80)
print("РЕЗУЛЬТАТЫ РАСЧЁТА В ТОЧКАХ ИЗМЕРЕНИЙ")
print("=" * 80)

print("\nОсновная таблица:")
print("f (МГц) | HPBW_H (\") | HPBW_V (\") | σ_H (\") | R_s/σ_H | Ω_A (сек²) | T_a (K)")
print("-" * 85)

for res in results:
    print(f"{res['freq']:7.0f} | {res['hpbw_h']:9.1f} | {res['hpbw_v']:9.0f} | "
          f"{res['sigma_h']:7.1f} | {res['R_s_sigma_h_ratio']:7.3f} | "
          f"{res['Omega_A']:9.2e} | {res['T_a']:7.1f}")

# ============================================================================
# 5. ДЕТАЛЬНЫЙ АНАЛИЗ ДЛЯ КРАЙНИХ ТОЧЕК
# ============================================================================
print("\n" + "=" * 80)
print("ДЕТАЛЬНЫЙ АНАЛИЗ ДЛЯ КРАЙНИХ ТОЧЕК")
print("=" * 80)

for f in [1000, 3000]:
    # Находим индекс ближайшей точки измерения
    idx = np.argmin(np.abs(freq_data - f))
    res = results[idx]

    print(f"\nЧастота {res['freq']:.0f} МГц:")
    print(f"  HPBW_H = {res['hpbw_h']:.1f}\"")
    print(f"  HPBW_V = {res['hpbw_v']:.0f}\"")
    print(f"  σ_H = {res['sigma_h']:.1f}\", σ_V = {res['sigma_v']:.0f}\"")
    print(f"  R_s/σ_H = {res['R_s_sigma_h_ratio']:.3f}")
    print(f"  I_src = {res['I_src']:.2e} угл.сек²")
    print(f"  Ω_A = {res['Omega_A']:.2e} угл.сек²")
    print(f"  I_src/Ω_A = {res['I_src'] / res['Omega_A']:.5f}")
    print(f"  T_a = {res['T_a']:.1f} K")

print(f"\nОтношение T_a({results[-1]['freq']:.0f})/T_a({results[0]['freq']:.0f}) = "
      f"{results[-1]['T_a'] / results[0]['T_a']:.2f}")

# ============================================================================
# 6. ГРАФИКИ (ТОЛЬКО ТОЧКИ ИЗМЕРЕНИЙ)
# ============================================================================
fig = plt.figure(figsize=(15, 10))

# График 1: HPBW_H (как на вашем графике - логарифмическая шкала по Y)
ax1 = plt.subplot(2, 3, 1)
ax1.scatter(freq_data, hpbw_h_data, color='orange', s=3, edgecolors='orange', linewidth=1.5)
ax1.set_xlabel('Частота, МГц', fontsize=12)
ax1.set_ylabel('HPBW_H, угл. секунд', fontsize=12)
ax1.set_title('Экспериментальные данные HPBW_H', fontsize=14)
ax1.set_yscale('log')
ax1.grid(True, alpha=0.3, which='both')
ax1.set_ylim(50, 300)

# Добавляем значения рядом с точками
# for i, (x, y) in enumerate(zip(freq_data, hpbw_h_data)):
#     ax1.annotate(f'{y:.0f}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=9)

# График 2: Антенная температура T_a
ax2 = plt.subplot(2, 3, 2)
T_a_values = [res['T_a'] for res in results]
ax2.scatter(freq_data, T_a_values, color='green', s=3, edgecolors='green', linewidth=1.5)
ax2.set_xlabel('Частота, МГц', fontsize=12)
ax2.set_ylabel('Антенная температура, K', fontsize=12)
ax2.set_title('Антенная температура в точках измерений', fontsize=14)
ax2.grid(True, alpha=0.3)

# Добавляем значения
# for i, (x, y) in enumerate(zip(freq_data, T_a_values)):
#     ax2.annotate(f'{y:.0f}', (x, y), xytext=(5, 5), textcoords='offset points', fontsize=9)

# График 3: Отношение R_s/σ_H
ax3 = plt.subplot(2, 3, 3)
ratios = [res['R_s_sigma_h_ratio'] for res in results]
ax3.scatter(freq_data, ratios, color='purple', s=3, edgecolors='black', linewidth=1.5)
ax3.axhline(y=1, color='red', linestyle='--', linewidth=1, label='R_s = σ_H')
ax3.set_xlabel('Частота, МГц', fontsize=12)
ax3.set_ylabel('R_s / σ_H', fontsize=12)
ax3.set_title('Отношение размера источника к ширине луча', fontsize=14)
ax3.grid(True, alpha=0.3)
ax3.legend()
ax3.set_ylim(0, 4)

# График 4: I_src и Ω_A
ax4 = plt.subplot(2, 3, 4)
I_src_vals = [res['I_src'] for res in results]
Omega_A_vals = [res['Omega_A'] for res in results]
ax4.scatter(freq_data, I_src_vals, color='blue', s=3, label='I_src', edgecolors='blue', linewidth=1.5)
ax4.scatter(freq_data, Omega_A_vals, color='red', s=3, label='Ω_A', edgecolors='red', linewidth=1.5)
ax4.set_xlabel('Частота, МГц', fontsize=12)
ax4.set_ylabel('Угловая площадь, сек²', fontsize=12)
ax4.set_title('Числитель (I_src) и знаменатель (Ω_A)', fontsize=14)
ax4.set_yscale('log')
ax4.grid(True, alpha=0.3, which='both')
ax4.legend()

# График 5: Отношение I_src/Ω_A
ax5 = plt.subplot(2, 3, 5)
ratios_io = [res['I_src'] / res['Omega_A'] for res in results]
ax5.scatter(freq_data, ratios_io, color='brown', s=3, edgecolors='black', linewidth=1.5)
ax5.set_xlabel('Частота, МГц', fontsize=12)
ax5.set_ylabel('I_src / Ω_A', fontsize=12)
ax5.set_title('Доля мощности источника в луче', fontsize=14)
ax5.grid(True, alpha=0.3)

# График 6: Сравнение HPBW_H и HPBW_V
ax6 = plt.subplot(2, 3, 6)
hpbw_v_vals = [res['hpbw_v'] for res in results]
ax6.scatter(freq_data, hpbw_h_data, color='orange', s=3, label='HPBW_H (эксп.)', edgecolors='black', linewidth=1.5)
ax6.scatter(freq_data, hpbw_v_vals, color='cyan', s=3, label='HPBW_V (∝ 1/f)', edgecolors='black', linewidth=1.5)
ax6.axhline(y=R_s_arcsec * 2, color='red', linestyle='--', linewidth=1, label='Диаметр источника')
ax6.set_xlabel('Частота, МГц', fontsize=12)
ax6.set_ylabel('HPBW, угл. секунд', fontsize=12)
ax6.set_title('Сравнение HPBW_H и HPBW_V', fontsize=14)
ax6.set_yscale('log')
ax6.grid(True, alpha=0.3, which='both')
ax6.legend()

plt.tight_layout()
plt.show()

# ============================================================================
# 7. ЭКСПОРТ РЕЗУЛЬТАТОВ В ФАЙЛ
# ============================================================================
import csv

# Сохранение результатов в CSV файл
output_filename = 'antenna_temperature_results.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ['freq_MHz', 'HPBW_H_arcsec', 'HPBW_V_arcsec',
                  'sigma_H_arcsec', 'sigma_V_arcsec', 'R_s_sigma_H_ratio',
                  'I_src_arcsec2', 'Omega_A_arcsec2', 'T_a_K']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for res in results:
        writer.writerow({
            'freq_MHz': res['freq'],
            'HPBW_H_arcsec': res['hpbw_h'],
            'HPBW_V_arcsec': res['hpbw_v'],
            'sigma_H_arcsec': res['sigma_h'],
            'sigma_V_arcsec': res['sigma_v'],
            'R_s_sigma_H_ratio': res['R_s_sigma_h_ratio'],
            'I_src_arcsec2': res['I_src'],
            'Omega_A_arcsec2': res['Omega_A'],
            'T_a_K': res['T_a']
        })

print(f"\nРезультаты сохранены в файл: {output_filename}")


# ============================================================================
# 8. ФУНКЦИЯ ДЛЯ РАСЧЁТА НА НОВЫХ ЧАСТОТАХ (БЕЗ ИНТЕРПОЛЯЦИИ!)
# ============================================================================
def calculate_at_measured_frequencies(new_freqs):
    """
    Расчёт антенной температуры только для тех частот,
    для которых есть экспериментальные данные HPBW_H
    """
    results_filtered = []

    for f in new_freqs:
        # Проверяем, есть ли эта частота в исходных данных
        if f in freq_data:
            idx = np.where(freq_data == f)[0][0]
            results_filtered.append(results[idx])
        else:
            print(f"Внимание: для частоты {f} МГц нет экспериментальных данных HPBW_H")

    return results_filtered


# Пример использования:
print("\n" + "=" * 80)
print("ПРИМЕР ИСПОЛЬЗОВАНИЯ ФУНКЦИИ ДЛЯ КОНКРЕТНЫХ ЧАСТОТ")
print("=" * 80)

# Запрашиваем T_a только для тех частот, где есть данные
requested_freqs = [1000, 1500, 2000, 2500, 3000]
filtered_results = calculate_at_measured_frequencies(requested_freqs)

for res in filtered_results:
    print(f"{res['freq']} МГц: T_a = {res['T_a']:.1f} K")

print("\n" + "=" * 80)
print("ВАЖНО: программа использует ТОЛЬКО экспериментальные точки,")
print("без интерполяции между ними!")
print("=" * 80)

# Как использовать:
# Замените данные в разделе 1 на ваши реальные измерения:
#
# python
# freq_data = ваш_массив_частот  # МГц
# hpbw_h_data = ваш_массив_HPBW_H  # угловых секунд
# Запустите программу — она рассчитает T_a только для этих частот.
#
# Если нужна T_a на частоте, которой нет в данных — добавьте эту частоту и соответствующее HPBW_H в массивы.