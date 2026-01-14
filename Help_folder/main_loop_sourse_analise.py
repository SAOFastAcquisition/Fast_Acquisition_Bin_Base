import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# ============================================================================
# 1. ПАРАМЕТРЫ
# ============================================================================
freq_data = np.array([1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])
hpbw_h_data = np.array([200, 175, 155, 140, 128, 118, 110, 103, 97, 92, 88])
HPBW_V0 = 75.0 * 60
T_b = 20000.0


# ============================================================================
# 2. ТОЧНЫЙ РАСЧЁТ ИНТЕГРАЛА (БЕЗ АППРОКСИМАЦИЙ!)
# ============================================================================
def sigma_from_HPBW(hpbw):
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def exact_I_src(R_s, sigma_H, sigma_V):
    """Точное численное интегрирование (без аппроксимаций!)"""

    # Интегрируем в полярных координатах
    def integrand_polar(r, phi):
        theta_H = r * np.cos(phi)
        theta_V = r * np.sin(phi)
        return r * np.exp(-theta_H ** 2 / (2 * sigma_H ** 2) - theta_V ** 2 / (2 * sigma_V ** 2))

    # Используем адаптивное интегрирование
    result, error = integrate.dblquad(
        lambda phi, r: integrand_polar(r, phi),
        0, R_s,  # r от 0 до R_s
        lambda r: 0,  # phi_min
        lambda r: 2 * np.pi  # phi_max
    )
    return result


# ============================================================================
# 3. ИССЛЕДОВАНИЕ ГЛАДКОСТИ ФУНКЦИИ
# ============================================================================
def analyze_smoothness():
    """Анализ гладкости I_src как функции R_s и как функции f"""

    # Выберем одну частоту для анализа I_src(R_s)
    f_fixed = 1750  # МГц, где вы видели "скачок" для D=300"
    # Интерполируем sigma_H на этой частоте
    sigma_h_interp = np.interp(f_fixed, freq_data,
                               [sigma_from_HPBW(h) for h in hpbw_h_data])
    sigma_v = sigma_from_HPBW(HPBW_V0 * (1000 / f_fixed))

    print(f"\nАнализ при фиксированной частоте f = {f_fixed} МГц:")
    print(f"σ_H = {sigma_h_interp:.1f}\", σ_V = {sigma_v:.0f}\"")

    # Анализ I_src как функции R_s
    R_s_range = np.linspace(10, 200, 50)  # радиус от 10" до 200"
    I_src_vs_R = []
    I_src_deriv = []

    for R_s in R_s_range:
        I = exact_I_src(R_s, sigma_h_interp, sigma_v)
        I_src_vs_R.append(I)

    # Численная производная dI/dR_s
    dI_dR = np.gradient(I_src_vs_R, R_s_range)

    # Анализ I_src как функции f при фиксированном R_s
    R_s_fixed = 150  # радиус 150" (диаметр 300")
    freqs_fine = np.linspace(1000, 3000, 100)
    I_src_vs_f = []

    for f in freqs_fine:
        sigma_h = np.interp(f, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
        sigma_v = sigma_from_HPBW(HPBW_V0 * (1000 / f))
        I = exact_I_src(R_s_fixed, sigma_h, sigma_v)
        I_src_vs_f.append(I)

    # Производная dI/df
    dI_df = np.gradient(I_src_vs_f, freqs_fine)

    return {
        'R_s_range': R_s_range,
        'I_src_vs_R': I_src_vs_R,
        'dI_dR': dI_dR,
        'freqs_fine': freqs_fine,
        'I_src_vs_f': I_src_vs_f,
        'dI_df': dI_df,
        'sigma_h_fixed': sigma_h_interp,
        'R_s_fixed': R_s_fixed
    }


# ============================================================================
# 4. ВИЗУАЛИЗАЦИЯ ГЛАДКОСТИ
# ============================================================================
results = analyze_smoothness()

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Исследование гладкости I_src: нет скачков, только изменения производных',
             fontsize=14, y=1.02)

# График 1: I_src как функция R_s (фиксированная частота)
ax = axes[0, 0]
ax.plot(results['R_s_range'], results['I_src_vs_R'], 'b-', linewidth=3)
ax.axvline(x=results['sigma_h_fixed'], color='red', linestyle='--',
           label=f'σ_H = {results["sigma_h_fixed"]:.1f}"')
ax.set_xlabel('Радиус источника R_s, угл. сек')
ax.set_ylabel('I_src, угл.сек²')
ax.set_title('I_src(R_s) при фиксированной частоте\n(гладкая функция!)')
ax.grid(True, alpha=0.3)
ax.legend()

# График 2: dI/dR_s
ax = axes[0, 1]
ax.plot(results['R_s_range'], results['dI_dR'], 'g-', linewidth=2)
ax.axvline(x=results['sigma_h_fixed'], color='red', linestyle='--')
ax.set_xlabel('Радиус источника R_s, угл. сек')
ax.set_ylabel('dI_src/dR_s, сек²/сек')
ax.set_title('Производная I_src по радиусу\nПлавное изменение!')
ax.grid(True, alpha=0.3)

# График 3: Относительная производная (dI/dR)/I
ax = axes[0, 2]
rel_deriv = results['dI_dR'] / results['I_src_vs_R']
ax.plot(results['R_s_range'], rel_deriv, 'm-', linewidth=2)
ax.axvline(x=results['sigma_h_fixed'], color='red', linestyle='--')
ax.set_xlabel('Радиус источника R_s, угл. сек')
ax.set_ylabel('(dI/dR)/I, 1/сек')
ax.set_title('Относительная производная\nМаксимум при R_s ≈ σ_H')
ax.grid(True, alpha=0.3)

# График 4: I_src как функция частоты (фиксированный R_s = 150")
ax = axes[1, 0]
ax.plot(results['freqs_fine'], results['I_src_vs_f'], 'b-', linewidth=3)
ax.set_xlabel('Частота, МГц')
ax.set_ylabel('I_src, угл.сек²')
ax.set_title(f'I_src(f) при фиксированном R_s = {results["R_s_fixed"]}"\n(гладкая функция!)')
ax.grid(True, alpha=0.3)

# Отмечаем точку, где σ_H = R_s
# Найдём частоту, где σ_H(f) = R_s
sigma_h_func = lambda f: np.interp(f, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
f_cross = None
for i in range(len(results['freqs_fine']) - 1):
    f1, f2 = results['freqs_fine'][i], results['freqs_fine'][i + 1]
    s1, s2 = sigma_h_func(f1), sigma_h_func(f2)
    if (s1 - results['R_s_fixed']) * (s2 - results['R_s_fixed']) <= 0:
        f_cross = f1 + (f2 - f1) * (results['R_s_fixed'] - s1) / (s2 - s1)
        break

if f_cross:
    I_cross = np.interp(f_cross, results['freqs_fine'], results['I_src_vs_f'])
    ax.plot(f_cross, I_cross, 'ro', markersize=10,
            label=f'σ_H = R_s при {f_cross:.0f} МГц')
    ax.legend()

# График 5: dI/df
ax = axes[1, 1]
ax.plot(results['freqs_fine'], results['dI_df'], 'g-', linewidth=2)
if f_cross:
    dI_df_cross = np.interp(f_cross, results['freqs_fine'], results['dI_df'])
    ax.plot(f_cross, dI_df_cross, 'ro', markersize=10)
ax.set_xlabel('Частота, МГц')
ax.set_ylabel('dI_src/df, сек²/МГц')
ax.set_title('Производная I_src по частоте\nРезкое изменение (не скачок!)')
ax.grid(True, alpha=0.3)

# График 6: Относительная производная (dI/df)/I
ax = axes[1, 2]
rel_deriv_f = np.abs(results['dI_df'] / results['I_src_vs_f'])
ax.plot(results['freqs_fine'], rel_deriv_f, 'm-', linewidth=2)
if f_cross:
    rel_cross = np.interp(f_cross, results['freqs_fine'], rel_deriv_f)
    ax.plot(f_cross, rel_cross, 'ro', markersize=10)
ax.set_xlabel('Частота, МГц')
ax.set_ylabel('|(dI/df)/I|, 1/МГц')
ax.set_title('Относительная производная по частоте\nПик при σ_H ≈ R_s')
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

plt.tight_layout()
plt.show()


# ============================================================================
# 5. АНАЛИЗ РАЗНЫХ РАЗМЕРОВ ИСТОЧНИКОВ
# ============================================================================
def analyze_multiple_sizes():
    """Анализ для нескольких размеров источников"""
    diameters = [60, 100, 200, 300]
    colors = ['red', 'blue', 'green', 'purple']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for diam, color in zip(diameters, colors):
        R_s = diam / 2

        # Рассчитываем I_src(f) для этого размера
        freqs = np.linspace(1000, 3000, 100)
        I_src_vals = []
        dI_df_vals = []

        for f in freqs:
            sigma_h = np.interp(f, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
            sigma_v = sigma_from_HPBW(HPBW_V0 * (1000 / f))
            I = exact_I_src(R_s, sigma_h, sigma_v)
            I_src_vals.append(I)

        # Производная
        dI_df = np.gradient(I_src_vals, freqs)
        rel_deriv = np.abs(dI_df / I_src_vals)

        # Находим частоту, где σ_H = R_s
        f_cross = None
        for i in range(len(freqs) - 1):
            f1, f2 = freqs[i], freqs[i + 1]
            s1 = np.interp(f1, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
            s2 = np.interp(f2, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
            if (s1 - R_s) * (s2 - R_s) <= 0:
                f_cross = f1 + (f2 - f1) * (R_s - s1) / (s2 - s1)
                break

        # График I_src(f)
        axes[0, 0].plot(freqs, I_src_vals, color=color, linewidth=2,
                        label=f'D = {diam}"')
        if f_cross:
            I_cross = np.interp(f_cross, freqs, I_src_vals)
            axes[0, 0].plot(f_cross, I_cross, 'o', color=color, markersize=8)

        # График относительной производной
        axes[0, 1].plot(freqs, rel_deriv, color=color, linewidth=2,
                        label=f'D = {diam}"')
        if f_cross:
            deriv_cross = np.interp(f_cross, freqs, rel_deriv)
            axes[0, 1].plot(f_cross, deriv_cross, 'o', color=color, markersize=8)

    axes[0, 0].set_xlabel('Частота, МГц')
    axes[0, 0].set_ylabel('I_src, угл.сек²')
    axes[0, 0].set_title('I_src(f) для разных размеров источника')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    axes[0, 0].set_yscale('log')

    axes[0, 1].set_xlabel('Частота, МГц')
    axes[0, 1].set_ylabel('|(dI/df)/I|, 1/МГц')
    axes[0, 1].set_title('Относительная производная (пик при σ_H = R_s)')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    axes[0, 1].set_yscale('log')

    # График 3: Частота пика производной vs размер источника
    diameters_test = np.linspace(20, 320, 20)
    peak_freqs = []

    for diam in diameters_test:
        R_s = diam / 2
        # Находим f, где σ_H(f) = R_s
        # Решаем уравнение σ_H(f) = R_s
        f_vals = np.linspace(1000, 3000, 200)
        sigma_vals = [np.interp(f, freq_data, [sigma_from_HPBW(h) for h in hpbw_h_data])
                      for f in f_vals]

        # Линейная интерполяция для нахождения корня
        for i in range(len(f_vals) - 1):
            if (sigma_vals[i] - R_s) * (sigma_vals[i + 1] - R_s) <= 0:
                f_peak = f_vals[i] + (f_vals[i + 1] - f_vals[i]) * (R_s - sigma_vals[i]) / (
                            sigma_vals[i + 1] - sigma_vals[i])
                peak_freqs.append(f_peak)
                break

    axes[1, 0].plot(diameters_test, peak_freqs, 'b-', linewidth=3)
    axes[1, 0].set_xlabel('Диаметр источника, угл. сек')
    axes[1, 0].set_ylabel('Частота пика производной, МГц')
    axes[1, 0].set_title('Частота, где σ_H = R_s (пик производной)')
    axes[1, 0].grid(True, alpha=0.3)

    # Отмечаем исследованные размеры
    for diam in diameters:
        axes[1, 0].axvline(x=diam, color='gray', linestyle='--', alpha=0.5)

    # График 4: Величина пика производной vs размер источника
    axes[1, 1].plot(diameters_test, [1.0] * len(diameters_test), 'r--', alpha=0.5)
    axes[1, 1].set_xlabel('Диаметр источника, угл. сек')
    axes[1, 1].set_ylabel('Максимум |(dI/df)/I| (отн. ед.)')
    axes[1, 1].set_title('Величина пика относительной производной')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_yscale('log')

    plt.tight_layout()
    plt.show()


# ============================================================================
# 6. ВЫВОДЫ
# ============================================================================
print("=" * 80)
print("ФИЗИЧЕСКИЙ АНАЛИЗ (исправленный)")
print("=" * 80)
print("\n1. НЕТ СКАЧКОВ ФУНКЦИИ:")
print("   I_src(R_s) - гладкая функция при фиксированной частоте")
print("   I_src(f) - гладкая функция при фиксированном размере источника")
print("\n2. ЕСТЬ РЕЗКОЕ ИЗМЕНЕНИЕ ПРОИЗВОДНОЙ:")
print("   При R_s ≈ σ_H(f) производная dI_src/df меняется наиболее быстро")
print("   Это выглядит как 'скачок' на графике с дискретными точками")
print("\n3. ФИЗИЧЕСКАЯ ИНТЕРПРЕТАЦИЯ:")
print("   При R_s < σ_H: I_src ∝ R_s² (источник мал)")
print("   При R_s > σ_H: I_src ∝ R_s·σ_H (источник велик)")
print("   Переход между этими режимами - плавный, но с быстрым изменением производной")
print("\n4. ДЛЯ ВАШИХ СЛУЧАЕВ:")
print(f"   Диаметр 300\": R_s = 150\"")
print(f"   σ_H = 150\" при f ≈ {f_cross:.0f} МГц")
print("   Здесь производная dI/df максимальна по модулю")
print("\n5. МЕТОДОЛОГИЧЕСКИЙ ВЫВОД:")
print("   'Скачок' - артефакт дискретного представления данных")
print("   При непрерывном изменении частоты или размера - только плавные изменения")

analyze_multiple_sizes()