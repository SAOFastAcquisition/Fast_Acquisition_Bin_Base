import numpy as np
import pickle
from pathlib import Path
import csv


#   Рабочая программа расчета антенной температуры, обусловленной активной областью
#
# ============================================================================
# 1. ЗАГРУЗКА ЭКСПЕРИМЕНТАЛЬНЫХ ДАННЫХ - ширина ДН по горизонтали на основе
#               измерений в автоколлимационном режиме
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
# Радиус активной области (диаметр 300")
R_s_arcsec = 100.0

# HPBW по вертикали на 1000 МГц (75 угловых минут)
HPBW_V0 = 75.0 * 60  # 4500 угл. секунд


# ============================================================================
# 3. ФОРМУЛА ЗИРИНА ДЛЯ ЯРКОСТНОЙ ТЕМПЕРАТУРЫ
# ============================================================================
def T_brightness_Zirin(f_MHz):
    """
    Яркостная температура активной области по Зирину
    T_b(f) = 140077 / (f/1000)^2.1 + 10880 [K]
    """
    return 140077.0 / (f_MHz / 1000.0) ** 2.1 + 10880.0


# ============================================================================
# 4. БЫСТРЫЙ РАСЧЁТ ИНТЕГРАЛА I_src (оптимизированная версия)
# ============================================================================
def sigma_from_HPBW(hpbw):
    """Преобразование HPBW в sigma для гауссова луча"""
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def I_src_fast(R_s, sigma_H, sigma_V, N_r=100, N_phi=200):
    """
    Быстрое численное интегрирование в полярных координатах
    Оптимизировано для скорости
    """
    # Предварительно вычисляем константы
    denom_H = 2.0 * sigma_H ** 2
    denom_V = 2.0 * sigma_V ** 2

    # Создаем сетки
    r = np.linspace(0, R_s, N_r)
    phi = np.linspace(0, 2 * np.pi, N_phi)

    dr = r[1] - r[0]
    dphi = phi[1] - phi[0]

    I_total = 0.0

    # Векторизованное вычисление (частично)
    for i in range(N_r - 1):
        ri = (r[i] + r[i + 1]) / 2.0

        # Используем векторизацию по углу
        phij = (phi[:-1] + phi[1:]) / 2.0

        # Декартовы координаты для всех углов одновременно
        theta_H = ri * np.cos(phij)
        theta_V = ri * np.sin(phij)

        # Гауссова ДН для всех углов
        P_n = np.exp(-theta_H ** 2 / denom_H - theta_V ** 2 / denom_V)

        # Суммируем вклад всех углов
        I_total += np.sum(P_n) * ri * dr * dphi

    return I_total


# ============================================================================
# 5. ОСНОВНАЯ ФУНКЦИЯ РАСЧЁТА
# ============================================================================
def calculate_antenna_temperature(freq_MHz, hpbw_h_arcsec, R_s_arcsec):
    """
    Расчёт антенной температуры для заданных параметров

    Возвращает: (T_b, T_a) где
    - T_b - яркостная температура по Зирину [K]
    - T_a - антенная температура [K]
    """
    # 1. Яркостная температура по Зирину
    T_b = T_brightness_Zirin(freq_MHz)

    # 2. HPBW по вертикали
    hpbw_v = HPBW_V0 * (1000.0 / freq_MHz)

    # 3. Стандартные отклонения гауссова луча
    sigma_h = sigma_from_HPBW(hpbw_h_arcsec)
    sigma_v = sigma_from_HPBW(hpbw_v)

    # 4. Интеграл по источнику
    I_src = I_src_fast(R_s_arcsec, sigma_h, sigma_v, N_r=100, N_phi=200)

    # 5. Площадь луча
    Omega_A = 2.0 * np.pi * sigma_h * sigma_v

    # 6. Антенная температура
    T_a = T_b * I_src / Omega_A

    return T_b, T_a


# ============================================================================
# 6. РАСЧЁТ ДЛЯ ВСЕХ ЧАСТОТ
# ============================================================================
print("\n" + "=" * 60)
print("РАСЧЁТ АНТЕННОЙ ТЕМПЕРАТУРЫ С T_b ПО ЗИРИНУ")
print("=" * 60)
print(f"Диаметр активной области: {2 * R_s_arcsec} угл. сек")
print(f"Количество точек: {len(freq_data)}")
print("\nЧастота (МГц) | T_b (K) | T_a (K) | T_a/T_b")
print("-" * 50)

results = []
for i in range(len(freq_data)):
    freq = freq_data[i]
    hpbw_h = hpbw_h_data[i]

    T_b, T_a = calculate_antenna_temperature(freq, hpbw_h, R_s_arcsec)
    efficiency = T_a / T_b

    results.append({
        'freq_MHz': freq,
        'T_b_K': T_b,
        'T_a_K': T_a,
        'capture_efficiency': efficiency,
        'HPBW_H_arcsec': hpbw_h
    })

    # Вывод только каждую 10-ю точку для краткости
    if i % max(1, len(freq_data) // 10) == 0 or i == len(freq_data) - 1:
        print(f"{freq:10.0f} | {T_b:7.0f} | {T_a:7.1f} | {efficiency:.4f}")

# ============================================================================
# 7. КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ
# ============================================================================
print("\n" + "=" * 60)
print("КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ")
print("=" * 60)

# Находим крайние частоты
min_freq = min(freq_data)
max_freq = max(freq_data)

# Находим соответствующие результаты
min_idx = np.argmin(np.abs([r['freq_MHz'] for r in results] - min_freq))
max_idx = np.argmin(np.abs([r['freq_MHz'] for r in results] - max_freq))

res_min = results[min_idx]
res_max = results[max_idx]

print(f"\nНа {min_freq:.0f} МГц:")
print(f"  T_b = {res_min['T_b_K']:.0f} K")
print(f"  T_a = {res_min['T_a_K']:.1f} K")
print(f"  Эффективность захвата = {res_min['capture_efficiency']:.4f}")

print(f"\nНа {max_freq:.0f} МГц:")
print(f"  T_b = {res_max['T_b_K']:.0f} K")
print(f"  T_a = {res_max['T_a_K']:.1f} K")
print(f"  Эффективность захвата = {res_max['capture_efficiency']:.4f}")

print(f"\nОтношения:")
print(f"  T_b({max_freq:.0f})/T_b({min_freq:.0f}) = {res_max['T_b_K'] / res_min['T_b_K']:.3f}")
print(f"  T_a({max_freq:.0f})/T_a({min_freq:.0f}) = {res_max['T_a_K'] / res_min['T_a_K']:.3f}")

# ============================================================================
# 8. ЭКСПОРТ В CSV
# ============================================================================
output_filename = f'antenna_temperature_Zirin_D{2 * R_s_arcsec:.0f}.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ['freq_MHz', 'HPBW_H_arcsec', 'T_b_K', 'T_a_K', 'capture_efficiency']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for res in results:
        writer.writerow({
            'freq_MHz': res['freq_MHz'],
            'HPBW_H_arcsec': res['HPBW_H_arcsec'],
            'T_b_K': res['T_b_K'],
            'T_a_K': res['T_a_K'],
            'capture_efficiency': res['capture_efficiency']
        })

print(f"\nРезультаты сохранены в файл: {output_filename}")


# ============================================================================
# 9. ОПЦИОНАЛЬНО: ПРОСТОЙ ГРАФИК
# ============================================================================
def plot_simple_results():
    """Простой график для визуализации"""
    import matplotlib.pyplot as plt

    freqs = [r['freq_MHz'] for r in results]
    T_b_vals = [r['T_b_K'] for r in results]
    T_a_vals = [r['T_a_K'] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # График 1: T_b и T_a
    ax1.scatter(freqs, T_b_vals, color='red', s=10, label='T_b по Зирину')
    ax1.scatter(freqs, T_a_vals, color='green', s=10, label='T_a')
    ax1.set_xlabel('Частота, МГц')
    ax1.set_ylabel('Температура, K')
    ax1.set_title(f'Яркостная и антенная температура\nДиаметр источника: {2 * R_s_arcsec}"')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_yscale('log')

    # График 2: Эффективность захвата
    efficiencies = [r['capture_efficiency'] for r in results]
    ax2.scatter(freqs, efficiencies, color='blue', s=10)
    ax2.set_xlabel('Частота, МГц')
    ax2.set_ylabel('T_a / T_b')
    ax2.set_title('Эффективность захвата источника')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# Спросить пользователя, нужен ли график
# plot_choice = input("\nПостроить простой график? (y/n): ").strip().lower()
# if plot_choice == 'y':
plot_simple_results()

print("\n" + "=" * 60)
print("РАСЧЁТ ЗАВЕРШЁН")
print("=" * 60)

# freq_MHz, HPBW_H_arcsec, T_b_K, T_a_K, capture_efficiency
# 1000, 200.0, 150000, 1234.5, 0.00823
# ...