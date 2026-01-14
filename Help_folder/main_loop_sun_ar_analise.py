import numpy as np
import pickle
from pathlib import Path
import csv
import matplotlib.pyplot as plt

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

print("Загружено точек:", len(freq_data))

# ============================================================================
# 2. ПАРАМЕТРЫ МОДЕЛИ
# ============================================================================
# Радиус активной области (диаметр 200")
R_ar_arcsec = 100.0  # 100" радиус = 200" диаметр

# Радиус Солнца (диаметр 960")
R_sun_arcsec = 480.0  # 480" радиус = 960" диаметр

# Положение AR относительно центра Солнца
AR_POS_X = 600.0  # угл. сек
AR_POS_Y = 0.0  # угл. сек

# Температурный фактор AR
T_AR_FACTOR = 10.0  # T_ar = 10 * T_quiet

# HPBW по вертикали на 1000 МГц (75 угловых минут)
HPBW_V0 = 75.0 * 60  # 4500 угл. секунд


# ============================================================================
# 3. ФОРМУЛА ЗИРИНА ДЛЯ ТЕМПЕРАТУРЫ СПОКОЙНОГО СОЛНЦА
# ============================================================================
def T_quiet_Sun_Zirin(f_MHz):
    """
    Яркостная температура спокойного Солнца по Зирину
    T_b(f) = 140077 / (f/1000)^2.1 + 10880 [K]
    """
    return 140077.0 / (f_MHz / 1000.0) ** 2.1 + 10880.0


# ============================================================================
# 4. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (ОСТАВЛЯЕМ КАК В СТАРОМ СКРИПТЕ)
# ============================================================================
def sigma_from_HPBW(hpbw):
    """Преобразование HPBW в sigma для гауссова луча"""
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def integrate_circle_fast(R, sigma_H, sigma_V, center_x=0.0, center_y=0.0, N_r=100, N_phi=200):
    """
    Быстрое интегрирование гауссова луча по кругу
    Аналогично I_src_fast из старого скрипта, но с поддержкой смещения центра

    Parameters:
    -----------
    R : float
        Радиус круга (источника)
    sigma_H, sigma_V : float
        Стандартные отклонения луча
    center_x, center_y : float
        Смещение центра круга относительно центра луча
    N_r, N_phi : int
        Количество точек интегрирования
    """
    # Предварительно вычисляем константы
    denom_H = 2.0 * sigma_H ** 2
    denom_V = 2.0 * sigma_V ** 2

    # Создаем сетки
    r = np.linspace(0, R, N_r)
    phi = np.linspace(0, 2 * np.pi, N_phi)

    dr = r[1] - r[0]
    dphi = phi[1] - phi[0]

    I_total = 0.0

    # Векторизованное вычисление
    for i in range(N_r - 1):
        ri = (r[i] + r[i + 1]) / 2.0

        # Используем векторизацию по углу
        phij = (phi[:-1] + phi[1:]) / 2.0

        # Декартовы координаты точек КРУГА
        x_circle = center_x + ri * np.cos(phij)
        y_circle = center_y + ri * np.sin(phij)

        # Гауссова ДН для всех углов
        P_n = np.exp(-x_circle ** 2 / denom_H - y_circle ** 2 / denom_V)

        # Суммируем вклад всех углов
        I_total += np.sum(P_n) * ri * dr * dphi

    return I_total


# ============================================================================
# 5. ФУНКЦИИ РАСЧЁТА (СОХРАНЯЕМ ЛОГИКУ СТАРОГО СКРИПТА)
# ============================================================================
def calculate_ar_only(freq_MHz, hpbw_h_arcsec, R_ar):
    """
    ТОЧНО ТАК ЖЕ, КАК В СТАРОМ СКРИПТЕ, но для AR в (0,0)
    Возвращает T_a только от AR (без фона)
    """
    # Температура AR (по Зирину)
    T_ar = T_quiet_Sun_Zirin(freq_MHz) * T_AR_FACTOR

    # HPBW по вертикали
    hpbw_v = HPBW_V0 * (1000.0 / freq_MHz)

    # Стандартные отклонения
    sigma_h = sigma_from_HPBW(hpbw_h_arcsec)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # Интеграл по AR (центр AR в (0,0) - как в старом скрипте)
    I_ar = integrate_circle_fast(R_ar, sigma_h, sigma_v, 0.0, 0.0)

    # Площадь луча
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # Антенная температура от AR
    T_a_ar = T_ar * I_ar / Omega_A

    return T_a_ar, I_ar, Omega_A


def calculate_with_background(freq_MHz, hpbw_h_arcsec, R_ar, R_sun, ar_x, ar_y):
    """
    Новая модель: Солнце + смещённая AR
    Центр луча всегда на AR (точка (ar_x, ar_y))
    """
    # Температуры
    T_quiet = T_quiet_Sun_Zirin(freq_MHz)
    T_ar = T_quiet * T_AR_FACTOR

    # HPBW по вертикали
    hpbw_v = HPBW_V0 * (1000.0 / freq_MHz)

    # Стандартные отклонения
    sigma_h = sigma_from_HPBW(hpbw_h_arcsec)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # Площадь луча
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # 1. Солнце: центр в (0,0), но луч смотрит на (ar_x, ar_y)
    #    => центр луча смещён относительно центра Солнца на (-ar_x, -ar_y)
    I_sun = integrate_circle_fast(R_sun, sigma_h, sigma_v, -ar_x, -ar_y)

    # 2. AR: центр AR в (ar_x, ar_y), луч смотрит туда же
    #    => центр AR совпадает с центром луча
    I_ar = integrate_circle_fast(R_ar, sigma_h, sigma_v, 0.0, 0.0)

    # Антенные температуры
    T_a_sun = T_quiet * I_sun / Omega_A
    T_a_ar = T_ar * I_ar / Omega_A
    T_a_total = T_a_sun + T_a_ar

    return {
        'T_quiet': T_quiet,
        'T_ar': T_ar,
        'T_a_sun': T_a_sun,
        'T_a_ar': T_a_ar,
        'T_a_total': T_a_total,
        'I_sun': I_sun,
        'I_ar': I_ar,
        'Omega_A': Omega_A,
        'sigma_h': sigma_h,
        'sigma_v': sigma_v
    }


# ============================================================================
# 6. ОСНОВНОЙ РАСЧЁТ И СРАВНЕНИЕ
# ============================================================================
print("\n" + "=" * 70)
print("СРАВНЕНИЕ: СТАРАЯ МОДЕЛЬ (AR в центре) vs НОВАЯ (Солнце + AR)")
print("=" * 70)
print(f"Диаметр AR: {2 * R_ar_arcsec} угл. сек")
print(f"Диаметр Солнца: {2 * R_sun_arcsec} угл. сек")
print(f"Положение AR: ({AR_POS_X}, {AR_POS_Y}) угл. сек")
print(f"T_ar = {T_AR_FACTOR} × T_quiet")

results_old = []  # Только AR в центре
results_new = []  # Солнце + смещённая AR

for i in range(len(freq_data)):
    freq = freq_data[i]
    hpbw_h = hpbw_h_data[i]

    # 1. Старая модель (только AR в центре)
    T_a_ar_old, I_ar_old, Omega_A_old = calculate_ar_only(freq, hpbw_h, R_ar_arcsec)

    # 2. Новая модель (Солнце + AR)
    res_new = calculate_with_background(freq, hpbw_h, R_ar_arcsec,
                                        R_sun_arcsec, AR_POS_X, AR_POS_Y)

    results_old.append({
        'freq': freq,
        'HPBW_H': hpbw_h,
        'T_a_ar': T_a_ar_old,
        'I_ar': I_ar_old,
        'Omega_A': Omega_A_old,
        'efficiency': I_ar_old / Omega_A_old
    })

    results_new.append({
        'freq': freq,
        'HPBW_H': hpbw_h,
        'T_quiet': res_new['T_quiet'],
        'T_ar': res_new['T_ar'],
        'T_a_sun': res_new['T_a_sun'],
        'T_a_ar': res_new['T_a_ar'],
        'T_a_total': res_new['T_a_total'],
        'I_sun': res_new['I_sun'],
        'I_ar': res_new['I_ar'],
        'Omega_A': res_new['Omega_A'],
        'sigma_h': res_new['sigma_h'],
        'sigma_v': res_new['sigma_v'],
        'ar_contribution': res_new['T_a_ar'] / res_new['T_a_total'] if res_new['T_a_total'] > 0 else 0
    })

    # Вывод прогресса
    if i % max(1, len(freq_data) // 20) == 0 or i == len(freq_data) - 1:
        print(f"{freq:6.0f} МГц | HPBW: {hpbw_h:5.1f}\" | "
              f"T_a_ar_old: {T_a_ar_old:7.1f} K | "
              f"T_a_ar_new: {res_new['T_a_ar']:7.1f} K | "
              f"T_a_total: {res_new['T_a_total']:7.1f} K")

# ============================================================================
# 7. АНАЛИЗ: ПОЧЕМУ T_a_ar_old И T_a_ar_new ДОЛЖНЫ БЫТЬ ОДИНАКОВЫМИ?
# ============================================================================
print("\n" + "=" * 70)
print("АНАЛИЗ РАСЧЁТА T_a_ar")
print("=" * 70)

print("\nВ старой модели (calculate_ar_only):")
print("  - Центр AR: (0, 0)")
print("  - Центр луча: (0, 0)")
print("  - I_ar = интеграл по кругу радиуса R_ar с центром в (0, 0)")

print(f"\nВ новой модели (calculate_with_background):")
print(f"  - Центр AR: ({AR_POS_X}, {AR_POS_Y})")
print(f"  - Центр луча: ({AR_POS_X}, {AR_POS_Y})")
print("  - I_ar = интеграл по кругу радиуса R_ar с центром в (0, 0) относительно луча")

print("\nВЕКТОР СМЕЩЕНИЯ:")
print(f"  Положение AR: ({AR_POS_X}, {AR_POS_Y})")
print(f"  Центр луча: ({AR_POS_X}, {AR_POS_Y})")
print(f"  => Относительное смещение: (0, 0)")

print("\nВЫВОД: T_a_ar_old и T_a_ar_new должны быть РАВНЫ!")
print("Потому что в обоих случаях:")
print("  1. Относительное положение AR и луча одинаково")
print("  2. Размеры AR одинаковы")
print("  3. Параметры луча одинаковы")

# Проверим это
print("\n" + "=" * 70)
print("ПРОВЕРКА РАВЕНСТВА T_a_ar_old и T_a_ar_new")
print("=" * 70)

max_diff = 0
max_diff_idx = 0
for i in range(len(results_old)):
    diff = abs(results_old[i]['T_a_ar'] - results_new[i]['T_a_ar'])
    if diff > max_diff:
        max_diff = diff
        max_diff_idx = i

print(f"Максимальная разница: {max_diff:.2e} K")
print(f"При частоте: {results_old[max_diff_idx]['freq']:.0f} МГц")
print(f"T_a_ar_old: {results_old[max_diff_idx]['T_a_ar']:.2f} K")
print(f"T_a_ar_new: {results_new[max_diff_idx]['T_a_ar']:.2f} K")

if max_diff < 1e-6:
    print("\n✓ T_a_ar_old и T_a_ar_new РАВНЫ с точностью до 1e-6 K")
else:
    print(f"\n⚠ Есть небольшая разница: {max_diff:.2e} K")
    print("  Возможные причины:")
    print("  1. Ошибки округления при интегрировании")
    print("  2. Разное количество точек интегрирования")

# ============================================================================
# 8. ГРАФИКИ ДЛЯ НАГЛЯДНОСТИ
# ============================================================================
fig = plt.figure(figsize=(15, 10))

# График 1: Сравнение T_a_ar
ax1 = plt.subplot(2, 3, 1)
freqs = [r['freq'] for r in results_old]
T_a_ar_old_vals = [r['T_a_ar'] for r in results_old]
T_a_ar_new_vals = [r['T_a_ar'] for r in results_new]

ax1.plot(freqs, T_a_ar_old_vals, 'b-', linewidth=2, label='T_a_ar (старая модель)')
ax1.plot(freqs, T_a_ar_new_vals, 'r--', linewidth=2, label='T_a_ar (новая модель)')
ax1.set_xlabel('Частота, МГц')
ax1.set_ylabel('T_a_ar, K')
ax1.set_title('Сравнение: T_a_ar от AR')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_yscale('log')

# График 2: Все компоненты новой модели
ax2 = plt.subplot(2, 3, 2)
T_a_sun_vals = [r['T_a_sun'] for r in results_new]
T_a_total_vals = [r['T_a_total'] for r in results_new]

ax2.plot(freqs, T_a_total_vals, 'k-', linewidth=2, label='T_a_total')
ax2.plot(freqs, T_a_sun_vals, 'orange', linewidth=2, label='T_a_sun (фон)')
ax2.plot(freqs, T_a_ar_new_vals, 'r-', linewidth=2, label='T_a_ar (AR)')
ax2.set_xlabel('Частота, МГц')
ax2.set_ylabel('Антенная температура, K')
ax2.set_title('Новая модель: все компоненты')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_yscale('log')

# График 3: Вклад AR в T_a_total
ax3 = plt.subplot(2, 3, 3)
ar_contrib = [r['ar_contribution'] * 100 for r in results_new]

ax3.plot(freqs, ar_contrib, 'purple', linewidth=2)
ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
ax3.set_xlabel('Частота, МГц')
ax3.set_ylabel('Вклад AR, %')
ax3.set_title(f'Вклад AR в T_a_total\n(AR в ({AR_POS_X},{AR_POS_Y}))')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(0, 100)

# График 4: Эффективность захвата
ax4 = plt.subplot(2, 3, 4)
eff_old = [r['efficiency'] for r in results_old]
eff_new = [results_new[i]['I_ar'] / results_new[i]['Omega_A'] for i in range(len(results_new))]

ax4.plot(freqs, eff_old, 'b-', linewidth=2, label='Только AR (старая)')
ax4.plot(freqs, eff_new, 'r--', linewidth=2, label='AR в новой модели')
ax4.set_xlabel('Частота, МГц')
ax4.set_ylabel('I_ar / Ω_A')
ax4.set_title('Эффективность захвата AR')
ax4.grid(True, alpha=0.3)
ax4.legend()
ax4.set_yscale('log')

# График 5: Ширина луча
ax5 = plt.subplot(2, 3, 5)
HPBW_H = [r['HPBW_H'] for r in results_old]
sigma_H = [r['sigma_h'] for r in results_new]
sigma_V = [r['sigma_v'] for r in results_new]

ax5.plot(freqs, HPBW_H, 'blue', linewidth=2, label='HPBW_H (измерение)')
ax5.plot(freqs, sigma_H, 'cyan', linewidth=1, label='σ_H')
ax5.plot(freqs, sigma_V, 'green', linewidth=1, label='σ_V')
ax5.axhline(y=R_ar_arcsec, color='red', linestyle='--',
            label=f'Радиус AR ({R_ar_arcsec}\")')
ax5.set_xlabel('Частота, МГц')
ax5.set_ylabel('Угловой размер, секунды')
ax5.set_title('Ширина луча')
ax5.grid(True, alpha=0.3)
ax5.legend()
ax5.set_yscale('log')

# График 6: Отношение T_a_ar_old / T_a_ar_new
ax6 = plt.subplot(2, 3, 6)
ratio = []
for i in range(len(results_old)):
    if results_new[i]['T_a_ar'] > 0:
        ratio.append(results_old[i]['T_a_ar'] / results_new[i]['T_a_ar'])
    else:
        ratio.append(1.0)

ax6.plot(freqs, ratio, 'black', linewidth=2)
ax6.axhline(y=1.0, color='red', linestyle='--', label='Идеальное равенство')
ax6.set_xlabel('Частота, МГц')
ax6.set_ylabel('T_a_ar_old / T_a_ar_new')
ax6.set_title('Отношение старый/новый расчет T_a_ar')
ax6.grid(True, alpha=0.3)
ax6.legend()
ax6.set_ylim(0.99, 1.01)

plt.tight_layout()
plt.show()

# ============================================================================
# 9. ЭКСПОРТ РЕЗУЛЬТАТОВ
# ============================================================================
output_filename = f'sun_AR_comparison_D{2 * R_ar_arcsec:.0f}.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = [
        'freq_MHz', 'HPBW_H_arcsec',
        'T_a_ar_old_K', 'T_a_ar_new_K', 'T_a_sun_K', 'T_a_total_K',
        'I_ar_old', 'I_ar_new', 'I_sun', 'Omega_A',
        'sigma_H', 'sigma_V', 'ar_contribution_%'
    ]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(len(results_old)):
        writer.writerow({
            'freq_MHz': results_old[i]['freq'],
            'HPBW_H_arcsec': results_old[i]['HPBW_H'],
            'T_a_ar_old_K': results_old[i]['T_a_ar'],
            'T_a_ar_new_K': results_new[i]['T_a_ar'],
            'T_a_sun_K': results_new[i]['T_a_sun'],
            'T_a_total_K': results_new[i]['T_a_total'],
            'I_ar_old': results_old[i]['I_ar'],
            'I_ar_new': results_new[i]['I_ar'],
            'I_sun': results_new[i]['I_sun'],
            'Omega_A': results_new[i]['Omega_A'],
            'sigma_H': results_new[i]['sigma_h'],
            'sigma_V': results_new[i]['sigma_v'],
            'ar_contribution_%': results_new[i]['ar_contribution'] * 100
        })

print(f"\nРезультаты сохранены в файл: {output_filename}")

# ============================================================================
# 10. ИТОГ
# ============================================================================
print("\n" + "=" * 70)
print("ИТОГОВЫЕ ВЫВОДЫ")
print("=" * 70)
print("\n1. T_a_ar (антенная температура от AR) НЕ МЕНЯЕТСЯ при добавлении фона.")
print("   Она зависит ТОЛЬКО от:")
print("   - Относительного положения AR и луча (в данном случае совпадают)")
print("   - Размера AR")
print("   - Параметров луча (HPBW_H, HPBW_V)")
print("   - Яркостной температуры AR")

print("\n2. Что меняется в новой модели:")
print("   - Добавляется T_a_sun (вклад фонового Солнца)")
print("   - Общая T_a_total = T_a_sun + T_a_ar")
print("   - Появляется понятие 'вклад AR в общий сигнал'")

print("\n3. Если в новой модели T_a_ar ведёт себя иначе:")
print("   ✓ Проверьте совпадение центров AR и луча")
print("   ✓ Убедитесь, что используется одинаковая формула для T_ar")
print("   ✓ Проверьте точность численного интегрирования")

print("\n" + "=" * 70)
print("РАСЧЁТ ЗАВЕРШЁН")
print("=" * 70)