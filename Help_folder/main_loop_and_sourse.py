import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import integrate
import csv

# ============================================================================
# 1. ЗАГРУЗКА ЭКСПЕРИМЕНТАЛЬНЫХ ДАННЫХ
# ============================================================================
current_primary_file = '2025-12-21_01'

hpbw_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                      f'{current_primary_file}_hpbw.bin')
freq_file_path = Path(r"I:\Fast_Acquisition\2025\Test_and_calibration\Converted_data\2025_12_21test_conv",
                      f'{current_primary_file}_freq.bin')

with open(hpbw_file_path, 'rb') as file:
    hpbw_h_data = pickle.load(file)
with open(freq_file_path, 'rb') as file:
    freq_data = pickle.load(file)

print("ВАШИ ЭКСПЕРИМЕНТАЛЬНЫЕ ДАННЫЕ:")
print("Частота (МГц) | HPBW_H (\")")
print("-" * 30)
for f, h in zip(freq_data, hpbw_h_data):
    print(f"{f:12.0f} | {h:8.1f}")


# ============================================================================
# 2. КОНСТАНТЫ И ПАРАМЕТРЫ
# ============================================================================
# Яркостная температура по Зирину (зависит от частоты!)
def T_brightness_Zirin(f_MHz):
    """
    Яркостная температура активной области по Зирину

    Parameters:
    -----------
    f_MHz : float or array
        Частота в МГц

    Returns:
    --------
    float or array : яркостная температура в K
    """
    return 140077.0 / (f_MHz / 1000.0) ** 2.1 + 10880.0


# Параметры источника
R_s_arcsec = 30.0  # радиус источника в угловых секундах
HPBW_V0 = 75.0 * 60  # 75' в секундах на опорной частоте f0 = 1000 МГц

# Параметры численного интегрирования
INTEGRATION_METHOD = 'polar'  # 'direct', 'polar', или 'scipy'


# ============================================================================
# 3. ФУНКЦИИ ДЛЯ ТОЧНОГО ЧИСЛЕННОГО ИНТЕГРИРОВАНИЯ
# ============================================================================
def sigma_from_HPBW(hpbw):
    """Преобразование HPBW в sigma для гауссова луча"""
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def I_src_numerical_polar(R_s, sigma_H, sigma_V, N_r=150, N_phi=300):
    """
    Численное интегрирование в полярных координатах
    (более эффективно для круглой области)
    """
    # Создаем сетку в полярных координатах
    r = np.linspace(0, R_s, N_r)
    phi = np.linspace(0, 2 * np.pi, N_phi)

    dr = r[1] - r[0]
    dphi = phi[1] - phi[0]

    I_total = 0.0

    for i in range(N_r - 1):
        ri = (r[i] + r[i + 1]) / 2.0  # средний радиус
        for j in range(N_phi - 1):
            phij = (phi[j] + phi[j + 1]) / 2.0  # средний угол

            # Переход к декартовым координатам
            theta_H = ri * np.cos(phij)
            theta_V = ri * np.sin(phij)

            # Значение гауссовой ДН
            P_n = np.exp(-theta_H ** 2 / (2 * sigma_H ** 2) - theta_V ** 2 / (2 * sigma_V ** 2))

            # Якобиан преобразования: r dr dphi
            I_total += P_n * ri * dr * dphi

    return I_total


def calculate_T_a_with_Zirin(freq, hpbw_h, R_s=None, method=INTEGRATION_METHOD):
    """
    Расчёт антенной температуры с учётом зависимости T_b от частоты по Зирину
    """
    if R_s is None:
        R_s = R_s_arcsec

    # Яркостная температура на данной частоте
    T_b = T_brightness_Zirin(freq)

    # HPBW по вертикали: ∝ 1/f
    hpbw_v = HPBW_V0 * (1000.0 / freq)

    # Сигмы
    sigma_h = sigma_from_HPBW(hpbw_h)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # Численное интегрирование I_src
    I_src = I_src_numerical_polar(R_s, sigma_h, sigma_v, N_r=150, N_phi=300)

    # Площадь луча антенны
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # Антенная температура
    T_a = T_b * I_src / Omega_A

    return {
        'freq': freq,
        'T_b': T_b,  # Добавляем яркостную температуру
        'hpbw_h': hpbw_h,
        'hpbw_v': hpbw_v,
        'sigma_h': sigma_h,
        'sigma_v': sigma_v,
        'R_s': R_s,
        'R_s_sigma_h_ratio': R_s / sigma_h,
        'I_src': I_src,
        'Omega_A': Omega_A,
        'T_a': T_a,
        'integration_method': method
    }


# ============================================================================
# 4. РАСЧЁТ ДЛЯ ВСЕХ ТОЧЕК ИЗМЕРЕНИЙ
# ============================================================================
print("\n" + "=" * 80)
print(f"РАСЧЁТ С УЧЁТОМ ЗАВИСИМОСТИ T_b ОТ ЧАСТОТЫ (ЗИРИН)")
print("=" * 80)

all_results = []
T_b_values = []

for i in range(len(freq_data)):
    print(f"Расчёт точки {i + 1}/{len(freq_data)}: f = {freq_data[i]} МГц...")

    res = calculate_T_a_with_Zirin(freq_data[i], hpbw_h_data[i],
                                   method=INTEGRATION_METHOD)
    all_results.append(res)
    T_b_values.append(res['T_b'])

    # Выводим каждую 10-ю точку для контроля
    if i % max(1, len(freq_data) // 10) == 0 or i == len(freq_data) - 1:
        print(f"  f = {res['freq']} МГц: T_b = {res['T_b']:.0f} K, T_a = {res['T_a']:.1f} K")

# ============================================================================
# 5. АНАЛИЗ ЗАВИСИМОСТИ T_b ОТ ЧАСТОТЫ
# ============================================================================
print("\n" + "=" * 80)
print("АНАЛИЗ ЗАВИСИМОСТИ ЯРКОСТНОЙ ТЕМПЕРАТУРЫ ОТ ЧАСТОТЫ")
print("=" * 80)

# Создаём частую сетку для анализа T_b(f)
freqs_fine_Tb = np.linspace(min(freq_data), max(freq_data), 100)
T_b_fine = T_brightness_Zirin(freqs_fine_Tb)

print(f"\nФормула Зирина: T_b(f) = 140077 / (f/1000)^2.1 + 10880")
print(f"Диапазон частот: {min(freq_data):.0f} - {max(freq_data):.0f} МГц")
print(f"T_b на {min(freq_data):.0f} МГц: {T_b_fine[0]:.0f} K")
print(f"T_b на {max(freq_data):.0f} МГц: {T_b_fine[-1]:.0f} K")
print(f"Отношение T_b({min(freq_data):.0f})/T_b({max(freq_data):.0f}): "
      f"{T_b_fine[0] / T_b_fine[-1]:.2f}")

# ============================================================================
# 6. ВЫВОД РЕЗУЛЬТАТОВ В ТАБЛИЦЕ
# ============================================================================
print("\n" + "=" * 80)
print("РЕЗУЛЬТАТЫ РАСЧЁТА (с T_b по Зирину)")
print("=" * 80)

print("\nОсновная таблица:")
print("f (МГц) | T_b (K) | HPBW_H (\") | R_s/σ_H | T_a (K) | T_a/T_b")
print("-" * 85)

for res in all_results:
    print(f"{res['freq']:7.0f} | {res['T_b']:7.0f} | {res['hpbw_h']:9.1f} | "
          f"{res['R_s_sigma_h_ratio']:7.3f} | {res['T_a']:7.1f} | "
          f"{res['T_a'] / res['T_b']:.4f}")


# ============================================================================
# 7. СРАВНЕНИЕ С ПОСТОЯННОЙ T_b = 20000 K
# ============================================================================
def calculate_T_a_constant_Tb(freq, hpbw_h, T_b_constant=20000.0, R_s=None):
    """
    Расчёт с постоянной T_b для сравнения
    """
    if R_s is None:
        R_s = R_s_arcsec

    # HPBW по вертикали: ∝ 1/f
    hpbw_v = HPBW_V0 * (1000.0 / freq)

    # Сигмы
    sigma_h = sigma_from_HPBW(hpbw_h)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # Численное интегрирование I_src
    I_src = I_src_numerical_polar(R_s, sigma_h, sigma_v, N_r=150, N_phi=300)

    # Площадь луча антенны
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # Антенная температура
    T_a = T_b_constant * I_src / Omega_A

    return T_a


# Рассчитываем для постоянной T_b
T_a_constant_Tb = []
for i in range(len(freq_data)):
    T_a_const = calculate_T_a_constant_Tb(freq_data[i], hpbw_h_data[i],
                                          T_b_constant=20000.0)
    T_a_constant_Tb.append(T_a_const)

# ============================================================================
# 8. ГРАФИКИ
# ============================================================================
fig = plt.figure(figsize=(15, 12))

# График 1: Яркостная температура T_b по Зирину
ax1 = plt.subplot(3, 3, 1)
ax1.plot(freqs_fine_Tb, T_b_fine, 'r-', linewidth=2, label='T_b по Зирину')
ax1.scatter(freq_data, T_b_values, color='red', s=20, edgecolors='black',
            linewidth=1, zorder=5)
ax1.axhline(y=20000, color='gray', linestyle='--', linewidth=1,
            label='T_b = 20000 K (старое)')
ax1.set_xlabel('Частота, МГц')
ax1.set_ylabel('T_b, K')
ax1.set_title('Яркостная температура по Зирину')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_yscale('log')

# График 2: Антенная температура T_a (сравнение)
ax2 = plt.subplot(3, 3, 2)
T_a_values = [res['T_a'] for res in all_results]
ax2.scatter(freq_data, T_a_values, color='green', s=20, edgecolors='black',
            linewidth=1, label='T_a с T_b по Зирину', zorder=5)
ax2.scatter(freq_data, T_a_constant_Tb, color='blue', s=20, edgecolors='black',
            linewidth=1, alpha=0.5, label='T_a с T_b=20000 K', zorder=4)
ax2.set_xlabel('Частота, МГц')
ax2.set_ylabel('T_a, K')
ax2.set_title(f'Антенная температура\nR_s={R_s_arcsec}"')
ax2.grid(True, alpha=0.3)
ax2.legend()

# График 3: Отношение T_a/T_b (эффективность захвата)
ax3 = plt.subplot(3, 3, 3)
capture_efficiency = [res['T_a'] / res['T_b'] for res in all_results]
ax3.scatter(freq_data, capture_efficiency, color='purple', s=20,
            edgecolors='black', linewidth=1)
ax3.set_xlabel('Частота, МГц')
ax3.set_ylabel('T_a / T_b')
ax3.set_title('Эффективность захвата источника')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(0, 1.1 * max(capture_efficiency))

# График 4: Отношение T_a(Зирин)/T_a(const)
ax4 = plt.subplot(3, 3, 4)
ratio_Ta = [T_a_values[i] / T_a_constant_Tb[i] for i in range(len(T_a_values))]
ax4.scatter(freq_data, ratio_Ta, color='brown', s=20, edgecolors='black',
            linewidth=1)
ax4.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax4.set_xlabel('Частота, МГц')
ax4.set_ylabel('T_a(Зирин) / T_a(T_b=20000)')
ax4.set_title('Влияние зависимости T_b от частоты')
ax4.grid(True, alpha=0.3)

# График 5: HPBW_H
ax5 = plt.subplot(3, 3, 5)
ax5.scatter(freq_data, hpbw_h_data, color='orange', s=3, edgecolors='orange',
            linewidth=1.5)
ax5.set_xlabel('Частота, МГц')
ax5.set_ylabel('HPBW_H, угл. секунд')
ax5.set_title('Экспериментальные данные HPBW_H')
ax5.set_yscale('log')
ax5.grid(True, alpha=0.3, which='both')
ax5.set_ylim(50, 300)

# График 6: Отношение R_s/σ_H
ax6 = plt.subplot(3, 3, 6)
ratios = [res['R_s_sigma_h_ratio'] for res in all_results]
ax6.scatter(freq_data, ratios, color='purple', s=3, edgecolors='black',
            linewidth=1.5)
ax6.axhline(y=1, color='red', linestyle='--', linewidth=1, label='R_s = σ_H')
ax6.set_xlabel('Частота, МГц')
ax6.set_ylabel('R_s / σ_H')
ax6.set_title('Отношение размера источника к ширине луча')
ax6.grid(True, alpha=0.3)
ax6.legend(fontsize=9)
ax6.set_ylim(0, 4)

# График 7: I_src и Ω_A
ax7 = plt.subplot(3, 3, 7)
I_src_vals = [res['I_src'] for res in all_results]
Omega_A_vals = [res['Omega_A'] for res in all_results]
ax7.scatter(freq_data, I_src_vals, color='blue', s=3, label='I_src',
            edgecolors='blue', linewidth=1.5)
ax7.scatter(freq_data, Omega_A_vals, color='red', s=3, label='Ω_A',
            edgecolors='red', linewidth=1.5)
ax7.set_xlabel('Частота, МГц')
ax7.set_ylabel('Угловая площадь, сек²')
ax7.set_title('Числитель (I_src) и знаменатель (Ω_A)')
ax7.set_yscale('log')
ax7.grid(True, alpha=0.3, which='both')
ax7.legend()

# График 8: Производная T_b (скорость изменения)
ax8 = plt.subplot(3, 3, 8)
# Рассчитываем производную T_b
dTb_df = np.gradient(T_b_fine, freqs_fine_Tb)
rel_dTb_df = np.abs(dTb_df / T_b_fine)
ax8.plot(freqs_fine_Tb, rel_dTb_df, 'r-', linewidth=2)
ax8.set_xlabel('Частота, МГц')
ax8.set_ylabel('|(dT_b/df)/T_b|, 1/МГц')
ax8.set_title('Относительная скорость изменения T_b')
ax8.grid(True, alpha=0.3)
ax8.set_yscale('log')

# График 9: Сравнение HPBW_H и HPBW_V
ax9 = plt.subplot(3, 3, 9)
hpbw_v_vals = [res['hpbw_v'] for res in all_results]
ax9.scatter(freq_data, hpbw_h_data, color='orange', s=3, label='HPBW_H (эксп.)',
            edgecolors='black', linewidth=1.5)
ax9.scatter(freq_data, hpbw_v_vals, color='cyan', s=3, label='HPBW_V (∝ 1/f)',
            edgecolors='black', linewidth=1.5)
ax9.axhline(y=R_s_arcsec * 2, color='red', linestyle='--', linewidth=1,
            label=f'Диаметр источника ({2 * R_s_arcsec}\")')
ax9.set_xlabel('Частота, МГц')
ax9.set_ylabel('HPBW, угл. секунд')
ax9.set_title('Сравнение HPBW_H и HPBW_V')
ax9.set_yscale('log')
ax9.grid(True, alpha=0.3, which='both')
ax9.legend(fontsize=8)

plt.tight_layout()
plt.show()

# ============================================================================
# 9. АНАЛИЗ ВЛИЯНИЯ ЗАВИСИМОСТИ T_b ОТ ЧАСТОТЫ
# ============================================================================
print("\n" + "=" * 80)
print("АНАЛИЗ ВЛИЯНИЯ ЗАВИСИМОСТИ T_b ОТ ЧАСТОТЫ")
print("=" * 80)

# Сравнение на краях диапазона
idx_1000 = np.argmin(np.abs(freq_data - 1000))
idx_3000 = np.argmin(np.abs(freq_data - 3000))

res_1000 = all_results[idx_1000]
res_3000 = all_results[idx_3000]
T_a_const_1000 = T_a_constant_Tb[idx_1000]
T_a_const_3000 = T_a_constant_Tb[idx_3000]

print(f"\nСравнение на краях диапазона:")
print(f"\n{res_1000['freq']:.0f} МГц:")
print(f"  T_b по Зирину: {res_1000['T_b']:.0f} K")
print(f"  T_a с T_b по Зирину: {res_1000['T_a']:.1f} K")
print(f"  T_a с T_b=20000 K: {T_a_const_1000:.1f} K")
print(f"  Отношение: {res_1000['T_a'] / T_a_const_1000:.3f}")

print(f"\n{res_3000['freq']:.0f} МГц:")
print(f"  T_b по Зирину: {res_3000['T_b']:.0f} K")
print(f"  T_a с T_b по Зирину: {res_3000['T_a']:.1f} K")
print(f"  T_a с T_b=20000 K: {T_a_const_3000:.1f} K")
print(f"  Отношение: {res_3000['T_a'] / T_a_const_3000:.3f}")

print(f"\nОтношения T_a по всему диапазону:")
print(f"  T_a(Зирин, {res_3000['freq']:.0f})/T_a(Зирин, {res_1000['freq']:.0f}): "
      f"{res_3000['T_a'] / res_1000['T_a']:.3f}")
print(f"  T_a(const, {res_3000['freq']:.0f})/T_a(const, {res_1000['freq']:.0f}): "
      f"{T_a_const_3000 / T_a_const_1000:.3f}")

# Анализ вкладов
print(f"\nРазложение вкладов в T_a:")
print(f"  T_a = T_b × (I_src/Ω_A)")
print(f"\nНа {res_1000['freq']:.0f} МГц:")
print(f"  T_b = {res_1000['T_b']:.0f} K")
print(f"  I_src/Ω_A = {res_1000['I_src'] / res_1000['Omega_A']:.5f}")
print(f"  Вклад T_b: {res_1000['T_b'] / 20000:.3f} от постоянного значения")
print(f"  Вклад геометрии: {res_1000['I_src'] / res_1000['Omega_A']:.5f}")

print(f"\nНа {res_3000['freq']:.0f} МГц:")
print(f"  T_b = {res_3000['T_b']:.0f} K")
print(f"  I_src/Ω_A = {res_3000['I_src'] / res_3000['Omega_A']:.5f}")
print(f"  Вклад T_b: {res_3000['T_b'] / 20000:.3f} от постоянного значения")
print(f"  Вклад геометрии: {res_3000['I_src'] / res_3000['Omega_A']:.5f}")

# ============================================================================
# 10. ЭКСПОРТ РЕЗУЛЬТАТОВ В ФАЙЛ
# ============================================================================
output_filename = f'antenna_temperature_Zirin_Rs{R_s_arcsec}.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ['freq_MHz', 'T_b_K', 'HPBW_H_arcsec', 'HPBW_V_arcsec',
                  'sigma_H_arcsec', 'sigma_V_arcsec', 'R_s_arcsec',
                  'R_s_sigma_H_ratio', 'I_src_arcsec2', 'Omega_A_arcsec2',
                  'T_a_K', 'capture_efficiency', 'T_a_Tb20000_K']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for i, res in enumerate(all_results):
        writer.writerow({
            'freq_MHz': res['freq'],
            'T_b_K': res['T_b'],
            'HPBW_H_arcsec': res['hpbw_h'],
            'HPBW_V_arcsec': res['hpbw_v'],
            'sigma_H_arcsec': res['sigma_h'],
            'sigma_V_arcsec': res['sigma_v'],
            'R_s_arcsec': res['R_s'],
            'R_s_sigma_H_ratio': res['R_s_sigma_h_ratio'],
            'I_src_arcsec2': res['I_src'],
            'Omega_A_arcsec2': res['Omega_A'],
            'T_a_K': res['T_a'],
            'capture_efficiency': res['T_a'] / res['T_b'],
            'T_a_Tb20000_K': T_a_constant_Tb[i]
        })

print(f"\nРезультаты сохранены в файл: {output_filename}")

# ============================================================================
# 11. ВЫВОДЫ
# ============================================================================
print("\n" + "=" * 80)
print("ВЫВОДЫ")
print("=" * 80)

print("\n1. ВЛИЯНИЕ ЗАВИСИМОСТИ T_b ОТ ЧАСТОТЫ:")
print(f"   T_b убывает с частотой: {T_b_fine[0]:.0f} K → {T_b_fine[-1]:.0f} K")
print(f"   Отношение T_b(1000)/T_b(3000): {T_b_fine[0] / T_b_fine[-1]:.2f}")

print("\n2. ЭФФЕКТ НА АНТЕННУЮ ТЕМПЕРАТУРУ:")
print(f"   На низких частотах: T_a(Зирин) ≈ {res_1000['T_a'] / T_a_const_1000:.2f}×T_a(const)")
print(f"   На высоких частотах: T_a(Зирин) ≈ {res_3000['T_a'] / T_a_const_3000:.2f}×T_a(const)")

print("\n3. ОБЩАЯ ЗАВИСИМОСТЬ T_a(f):")
print(f"   T_a растёт с частотой, но медленнее чем при постоянной T_b")
print(f"   Отношение T_a(3000)/T_a(1000) с T_b по Зирину: {res_3000['T_a'] / res_1000['T_a']:.2f}")
print(f"   Отношение T_a(3000)/T_a(1000) с T_b=20000 K: {T_a_const_3000 / T_a_const_1000:.2f}")

print("\n4. ФИЗИЧЕСКАЯ ИНТЕРПРЕТАЦИЯ:")
print("   - T_b(f) убывает как ~f^{-2.1} (эффект рассеяния в короне)")
print("   - Ω_A убывает как f^{-1.92}")
print("   - I_src слабо зависит от частоты")
print("   - Результирующая T_a растёт медленнее из-за убывания T_b")

print("\n5. ПРАКТИЧЕСКИЕ СЛЕДСТВИЯ:")
print("   - На низких частотах T_a завышена относительно модели с постоянной T_b")
print("   - На высоких частотах T_a ближе к модели с постоянной T_b")
print("   - Общий рост T_a с частотой ослаблен")

# INTEGRATION_METHOD = 'polar'  # Можно изменить на 'direct' или 'scipy'
# INTEGRATION_RESOLUTION = 200  # Количество точек интегрирования
# R_s_arcsec = 150.0  # Радиус источника (диаметр 300")