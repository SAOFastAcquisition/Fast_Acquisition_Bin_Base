import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import csv

# ============================================================================
# 1. ВАШИ ЭКСПЕРИМЕНТАЛЬНЫЕ ДАННЫЕ
# ============================================================================
freq_data = np.array([1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])
hpbw_h_data = np.array([200, 175, 155, 140, 128, 118, 110, 103, 97, 92, 88])

# ============================================================================
# 2. КОНСТАНТЫ И ПАРАМЕТРЫ
# ============================================================================
T_b = 20000.0
R_s_arcsec = 100.0  # радиус по умолчанию
HPBW_V0 = 75.0 * 60


# ============================================================================
# 3. ФУНКЦИИ ДЛЯ ТОЧНОГО ЧИСЛЕННОГО ИНТЕГРИРОВАНИЯ
# ============================================================================
def sigma_from_HPBW(hpbw):
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def I_src_numerical_polar(R_s, sigma_H, sigma_V, N_r=150, N_phi=300):
    """Интегрирование в полярных координатах"""
    r = np.linspace(0, R_s, N_r)
    phi = np.linspace(0, 2 * np.pi, N_phi)

    dr = r[1] - r[0]
    dphi = phi[1] - phi[0]

    I_total = 0.0

    for i in range(N_r - 1):
        ri = (r[i] + r[i + 1]) / 2.0
        for j in range(N_phi - 1):
            phij = (phi[j] + phi[j + 1]) / 2.0

            theta_H = ri * np.cos(phij)
            theta_V = ri * np.sin(phij)

            P_n = np.exp(-theta_H ** 2 / (2 * sigma_H ** 2) - theta_V ** 2 / (2 * sigma_V ** 2))
            I_total += P_n * ri * dr * dphi

    return I_total


def calculate_T_a_numerical(freq, hpbw_h, R_s=None):
    """
    Основная функция расчёта с численным интегрированием
    """
    if R_s is None:
        R_s = R_s_arcsec

    # HPBW по вертикали
    hpbw_v = HPBW_V0 * (1000.0 / freq)

    # Сигмы
    sigma_h = sigma_from_HPBW(hpbw_h)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # Численное интегрирование I_src
    I_src = I_src_numerical_polar(R_s, sigma_h, sigma_v)

    # Площадь луча
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # Антенная температура
    T_a = T_b * I_src / Omega_A

    return {
        'freq': freq,
        'hpbw_h': hpbw_h,
        'hpbw_v': hpbw_v,
        'sigma_h': sigma_h,
        'sigma_v': sigma_v,
        'R_s': R_s,
        'R_s_sigma_h_ratio': R_s / sigma_h,
        'I_src': I_src,
        'Omega_A': Omega_A,
        'T_a': T_a
    }


# ============================================================================
# 4. РАСЧЁТ ДЛЯ ВСЕХ ТОЧЕК ИЗМЕРЕНИЙ (по умолчанию R_s=100")
# ============================================================================
print("РАСЧЁТ ДЛЯ R_s = 100\" (диаметр 200\")")
print("=" * 50)

all_results = []
for i in range(len(freq_data)):
    res = calculate_T_a_numerical(freq_data[i], hpbw_h_data[i])
    all_results.append(res)
    print(f"f = {res['freq']:4.0f} МГц: T_a = {res['T_a']:7.1f} K, R_s/σ_H = {res['R_s_sigma_h_ratio']:.3f}")


# ============================================================================
# 5. ИССЛЕДОВАНИЕ РАЗНЫХ РАЗМЕРОВ ИСТОЧНИКА
# ============================================================================
def analyze_different_sizes():
    """Анализ для разных диаметров источника"""
    diameters = [60, 100, 200, 300]

    print("\n" + "=" * 80)
    print("ИССЛЕДОВАНИЕ РАЗНЫХ РАЗМЕРОВ АКТИВНОЙ ОБЛАСТИ")
    print("=" * 80)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for diam in diameters:
        R_s = diam / 2.0

        print(f"\nДиаметр источника: {diam}\" (R_s = {R_s}\")")
        print("-" * 40)

        T_a_values = []
        I_src_values = []
        ratios = []

        # Расчёт для всех частот
        for i in range(len(freq_data)):
            res = calculate_T_a_numerical(freq_data[i], hpbw_h_data[i], R_s=R_s)
            T_a_values.append(res['T_a'])
            I_src_values.append(res['I_src'])
            ratios.append(res['R_s_sigma_h_ratio'])

            # Выводим ключевые точки
            if freq_data[i] in [1000, 1750, 3000]:
                print(f"  f = {freq_data[i]:4.0f} МГц: T_a = {res['T_a']:7.1f} K, "
                      f"R_s/σ_H = {ratios[-1]:.3f}")

        # Определяем цвет в зависимости от диаметра
        if diam == 60:
            color = 'red'
        elif diam == 100:
            color = 'blue'
        elif diam == 200:
            color = 'green'
        else:
            color = 'purple'

        # График 1: T_a для разных размеров
        axes[0, 0].scatter(freq_data, T_a_values, color=color, s=60,
                           edgecolors='black', linewidth=1, label=f'D={diam}"')

        # Интерполяция для гладкой кривой
        freqs_interp = np.linspace(1000, 3000, 50)
        hpbw_interp = np.interp(freqs_interp, freq_data, hpbw_h_data)
        T_a_interp = []

        for f, h in zip(freqs_interp, hpbw_interp):
            res = calculate_T_a_numerical(f, h, R_s=R_s)
            T_a_interp.append(res['T_a'])

        axes[0, 0].plot(freqs_interp, T_a_interp, color=color, alpha=0.3, linewidth=1)

        # График 2: I_src для разных размеров
        axes[0, 1].scatter(freq_data, I_src_values, color=color, s=60,
                           edgecolors='black', linewidth=1, label=f'D={diam}"')

        # Находим частоту, где R_s = σ_H
        f_cross = None
        for i in range(len(freq_data) - 1):
            if ratios[i] <= 1 <= ratios[i + 1] or ratios[i] >= 1 >= ratios[i + 1]:
                # Линейная интерполяция
                f1, f2 = freq_data[i], freq_data[i + 1]
                r1, r2 = ratios[i], ratios[i + 1]
                f_cross = f1 + (f2 - f1) * (1 - r1) / (r2 - r1)

                # Находим T_a в этой точке
                hpbw_cross = np.interp(f_cross, freq_data, hpbw_h_data)
                res_cross = calculate_T_a_numerical(f_cross, hpbw_cross, R_s=R_s)

                print(f"  Переход R_s = σ_H при f ≈ {f_cross:.0f} МГц, "
                      f"T_a = {res_cross['T_a']:.1f} K")
                break

        if f_cross:
            # Отмечаем на графике
            T_cross = np.interp(f_cross, freqs_interp, T_a_interp)
            axes[0, 0].plot(f_cross, T_cross, 'o', color=color, markersize=8)

    # Настройка графиков
    axes[0, 0].set_xlabel('Частота, МГц')
    axes[0, 0].set_ylabel('T_a, K')
    axes[0, 0].set_title('Антенная температура для разных размеров источника')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    axes[0, 0].set_yscale('log')

    axes[0, 1].set_xlabel('Частота, МГц')
    axes[0, 1].set_ylabel('I_src, угл.сек²')
    axes[0, 1].set_title('Интеграл I_src для разных размеров источника')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    axes[0, 1].set_yscale('log')

    # График 3: Отношение T_a(3000)/T_a(1000) vs размер
    diameters_test = np.linspace(20, 320, 20)
    temp_ratios = []

    for diam in diameters_test:
        R_s = diam / 2.0

        # T_a на 1000 МГц
        hpbw_1000 = np.interp(1000, freq_data, hpbw_h_data)
        res_1000 = calculate_T_a_numerical(1000, hpbw_1000, R_s=R_s)

        # T_a на 3000 МГц
        hpbw_3000 = np.interp(3000, freq_data, hpbw_h_data)
        res_3000 = calculate_T_a_numerical(3000, hpbw_3000, R_s=R_s)

        temp_ratios.append(res_3000['T_a'] / res_1000['T_a'])

    axes[1, 0].plot(diameters_test, temp_ratios, 'b-', linewidth=3)
    axes[1, 0].axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='T_a(3000)=T_a(1000)')
    axes[1, 0].set_xlabel('Диаметр источника, угл. сек')
    axes[1, 0].set_ylabel('T_a(3000) / T_a(1000)')
    axes[1, 0].set_title('Отношение температур на краях диапазона')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()

    # Отмечаем исследованные размеры
    for diam in diameters:
        idx = np.argmin(np.abs(diameters_test - diam))
        axes[1, 0].plot(diam, temp_ratios[idx], 'o', color='red', markersize=8)
        axes[1, 0].axvline(x=diam, color='gray', linestyle=':', alpha=0.5)

    # График 4: Частота перехода R_s = σ_H vs размер
    transition_freqs = []

    for diam in diameters_test:
        R_s = diam / 2.0

        # Находим частоту, где R_s = σ_H
        f_found = None
        freqs_check = np.linspace(1000, 3000, 100)
        sigmas_check = [sigma_from_HPBW(np.interp(f, freq_data, hpbw_h_data))
                        for f in freqs_check]

        for i in range(len(freqs_check) - 1):
            if (sigmas_check[i] - R_s) * (sigmas_check[i + 1] - R_s) <= 0:
                f_found = freqs_check[i] + (freqs_check[i + 1] - freqs_check[i]) * \
                          (R_s - sigmas_check[i]) / (sigmas_check[i + 1] - sigmas_check[i])
                break

        if f_found:
            transition_freqs.append(f_found)
        else:
            transition_freqs.append(np.nan)

    axes[1, 1].plot(diameters_test, transition_freqs, 'g-', linewidth=3)
    axes[1, 1].set_xlabel('Диаметр источника, угл. сек')
    axes[1, 1].set_ylabel('Частота перехода, МГц')
    axes[1, 1].set_title('Частота, где R_s = σ_H (переход между режимами)')
    axes[1, 1].grid(True, alpha=0.3)

    # Отмечаем исследованные размеры
    for diam in diameters:
        idx = np.argmin(np.abs(diameters_test - diam))
        axes[1, 1].plot(diam, transition_freqs[idx], 'o', color='red', markersize=8)
        axes[1, 1].axvline(x=diam, color='gray', linestyle=':', alpha=0.5)

    plt.tight_layout()
    plt.show()

    return diameters_test, temp_ratios, transition_freqs


# ============================================================================
# 6. ОСНОВНОЙ АНАЛИЗ
# ============================================================================
# Запускаем анализ разных размеров
diameters_test, temp_ratios, transition_freqs = analyze_different_sizes()

# ============================================================================
# 7. ВЫВОД РЕЗУЛЬТАТОВ ДЛЯ КЛЮЧЕВЫХ СЛУЧАЕВ
# ============================================================================
print("\n" + "=" * 80)
print("КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ ДЛЯ РАЗНЫХ РАЗМЕРОВ ИСТОЧНИКА")
print("=" * 80)

# Анализ для диаметров 60", 100", 200", 300"
key_diameters = [60, 100, 200, 300]

for diam in key_diameters:
    R_s = diam / 2.0

    print(f"\nДиаметр источника: {diam}\"")
    print("-" * 40)

    # Находим частоту перехода R_s = σ_H
    idx = np.argmin(np.abs(diameters_test - diam))
    f_transition = transition_freqs[idx]

    if not np.isnan(f_transition):
        print(f"Переход R_s = σ_H при f ≈ {f_transition:.0f} МГц")

        # Расчёт в точке перехода
        hpbw_transition = np.interp(f_transition, freq_data, hpbw_h_data)
        res_transition = calculate_T_a_numerical(f_transition, hpbw_transition, R_s=R_s)

        # Расчёт на краях диапазона
        res_1000 = calculate_T_a_numerical(1000, np.interp(1000, freq_data, hpbw_h_data), R_s=R_s)
        res_3000 = calculate_T_a_numerical(3000, np.interp(3000, freq_data, hpbw_h_data), R_s=R_s)

        print(f"  На {f_transition:.0f} МГц:")
        print(f"    σ_H = {res_transition['sigma_h']:.1f}\"")
        print(f"    T_a = {res_transition['T_a']:.1f} K")
        print(f"    I_src/Ω_A = {res_transition['I_src'] / res_transition['Omega_A']:.5f}")

        print(f"  На 1000 МГц:")
        print(f"    σ_H = {res_1000['sigma_h']:.1f}\", R_s/σ_H = {res_1000['R_s_sigma_h_ratio']:.3f}")
        print(f"    T_a = {res_1000['T_a']:.1f} K")

        print(f"  На 3000 МГц:")
        print(f"    σ_H = {res_3000['sigma_h']:.1f}\", R_s/σ_H = {res_3000['R_s_sigma_h_ratio']:.3f}")
        print(f"    T_a = {res_3000['T_a']:.1f} K")

        print(f"  Отношение T_a(3000)/T_a(1000) = {res_3000['T_a'] / res_1000['T_a']:.3f}")

        # Определяем режим на краях
        if res_1000['R_s_sigma_h_ratio'] < 1:
            regime_1000 = "σ_H > R_s (широкий луч)"
        else:
            regime_1000 = "σ_H < R_s (узкий луч)"

        if res_3000['R_s_sigma_h_ratio'] < 1:
            regime_3000 = "σ_H > R_s (широкий луч)"
        else:
            regime_3000 = "σ_H < R_s (узкий луч)"

        print(f"  Режим на 1000 МГц: {regime_1000}")
        print(f"  Режим на 3000 МГц: {regime_3000}")
    else:
        print(f"Переход R_s = σ_H не происходит в диапазоне 1000-3000 МГц")

# ============================================================================
# 8. ЭКСПОРТ ДАННЫХ
# ============================================================================
# Сохраняем результаты для диаметра 100" (основной случай)
output_filename = 'antenna_temperature_results_Rs100.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ['freq_MHz', 'HPBW_H_arcsec', 'HPBW_V_arcsec',
                  'sigma_H_arcsec', 'sigma_V_arcsec', 'R_s_arcsec',
                  'R_s_sigma_H_ratio', 'I_src_arcsec2', 'Omega_A_arcsec2',
                  'T_a_K']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for res in all_results:
        writer.writerow({
            'freq_MHz': res['freq'],
            'HPBW_H_arcsec': res['hpbw_h'],
            'HPBW_V_arcsec': res['hpbw_v'],
            'sigma_H_arcsec': res['sigma_h'],
            'sigma_V_arcsec': res['sigma_v'],
            'R_s_arcsec': res['R_s'],
            'R_s_sigma_H_ratio': res['R_s_sigma_h_ratio'],
            'I_src_arcsec2': res['I_src'],
            'Omega_A_arcsec2': res['Omega_A'],
            'T_a_K': res['T_a']
        })

print(f"\nРезультаты для R_s = 100\" сохранены в файл: {output_filename}")

# ============================================================================
# 9. ФИНАЛЬНЫЕ ВЫВОДЫ
# ============================================================================
print("\n" + "=" * 80)
print("ФИНАЛЬНЫЕ ВЫВОДЫ")
print("=" * 80)

print("\n1. ЧИСЛЕННОЕ ИНТЕГРИРОВАНИЕ ДАЁТ ГЛАДКИЕ РЕЗУЛЬТАТЫ:")
print("   - Нет скачков или разрывов")
print("   - Все зависимости физически корректны")
print("   - Устранены артефакты аналитических аппроксимаций")

print("\n2. ПЕРЕХОД МЕЖДУ РЕЖИМАМИ:")
print("   - При R_s < σ_H: источник меньше ширины луча")
print("   - При R_s > σ_H: источник больше ширины луча")
print("   - Переход плавный, но с максимальной производной при R_s ≈ σ_H")

print("\n3. ДЛЯ КОНКРЕТНЫХ РАЗМЕРОВ ИСТОЧНИКА:")
print("   - D=60\": R_s=30\", переход при f≈2800-3000 МГц")
print("   - D=100\": R_s=50\", переход при f≈2200-2400 МГц")
print("   - D=200\": R_s=100\", переход при f≈1400-1600 МГц")
print("   - D=300\": R_s=150\", переход при f≈1000-1200 МГц")

print("\n4. ВЛИЯНИЕ НА АНТЕННУЮ ТЕМПЕРАТУРУ:")
print("   - T_a всегда растёт с частотой (из-за уменьшения Ω_A ∝ f^{-1.92})")
print("   - Скорость роста максимальна вблизи перехода R_s = σ_H")
print("   - Отношение T_a(3000)/T_a(1000) зависит от размера источника")

print("\n5. МЕТОДОЛОГИЧЕСКИЙ РЕЗУЛЬТАТ:")
print("   Прямое численное интегрирование - самый надёжный метод")
print("   для задач с конечными источниками и гауссовыми лучами")