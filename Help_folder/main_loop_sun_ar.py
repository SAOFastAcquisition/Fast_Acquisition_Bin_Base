import numpy as np
import pickle
from pathlib import Path
import csv
import matplotlib.pyplot as plt

#   Рабочая программа расчета антенной температуры, обусловленной активной областью и спокойным Солнцем
#
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
R_sun_arcsec = 960.0  # 480" радиус = 960" диаметр

# Положение AR относительно центра Солнца
AR_POS_X = 300.0  # угл. сек
AR_POS_Y = 0.0  # угл. сек

# Температурный фактор AR
T_AR_FACTOR = 100.0  # T_ar = 10 * T_quiet

# HPBW по вертикали на 1000 МГц (75 угловых минут)
HPBW_V0 = 215.0 * 60  # 4500 угл. секунд


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
# 4. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================================
def sigma_from_HPBW(hpbw):
    """Преобразование HPBW в sigma для гауссова луча"""
    return hpbw / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def integrate_circle_fast(R, sigma_H, sigma_V, center_x=0.0, center_y=0.0, N_r=100, N_phi=200):
    """
    Быстрое интегрирование гауссова луча по кругу
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
# 5. ФУНКЦИЯ РАСЧЁТА С ФОНОМ (ИСПРАВЛЕННАЯ)
# ============================================================================
def calculate_with_background(freq_MHz, hpbw_h_arcsec, R_ar, R_sun, ar_x, ar_y):
    """
    Новая модель: Солнце + смещённая AR
    Центр луча ВСЕГДА в (ar_x, 0) - горизонтальная координата как у AR, но θ_y = 0
    AR может иметь произвольные координаты (ar_x, ar_y)
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

    # 1. Солнце: центр в (0, 0), луч в (ar_x, 0)
    #    => относительное смещение (-ar_x, 0)
    I_sun = integrate_circle_fast(R_sun, sigma_h, sigma_v, -ar_x, 0.0)

    # 2. AR: центр AR в (ar_x, ar_y), луч в (ar_x, 0)
    #    => относительное смещение (0, -ar_y)
    I_ar = integrate_circle_fast(R_ar, sigma_h, sigma_v, 0.0, -ar_y)

    # Антенные температуры
    T_a_sun = T_quiet * I_sun / Omega_A
    T_a_ar = T_ar * I_ar / Omega_A
    T_a_total = T_a_sun + T_a_ar

    return {
        'freq_MHz': freq_MHz,
        'HPBW_H': hpbw_h_arcsec,
        'T_quiet': T_quiet,
        'T_ar': T_ar,
        'T_a_sun': T_a_sun,
        'T_a_ar': T_a_ar,
        'T_a_total': T_a_total,
        'I_sun': I_sun,
        'I_ar': I_ar,
        'Omega_A': Omega_A,
        'sigma_h': sigma_h,
        'sigma_v': sigma_v,
        'ar_contribution': T_a_ar / T_a_total if T_a_total > 0 else 0
    }

# ============================================================================
# 6. РАСЧЁТ ДЛЯ ВСЕХ ЧАСТОТ
# ============================================================================
print("\n" + "=" * 70)
print("РАСЧЁТ МОДЕЛИ: СПОКОЙНОЕ СОЛНЦЕ + АКТИВНАЯ ОБЛАСТЬ")
print("=" * 70)
print(f"Диаметр AR: {2 * R_ar_arcsec} угл. сек")
print(f"Диаметр Солнца: {2 * R_sun_arcsec} угл. сек")
print(f"Положение AR: ({AR_POS_X}, {AR_POS_Y}) угл. сек")
print(f"T_ar = {T_AR_FACTOR} × T_quiet")

results = []

for i in range(len(freq_data)):
    freq = freq_data[i]
    hpbw_h = hpbw_h_data[i]

    res = calculate_with_background(freq, hpbw_h, R_ar_arcsec,
                                    R_sun_arcsec, AR_POS_X, AR_POS_Y)
    results.append(res)

    # Вывод прогресса
    if i % max(1, len(freq_data) // 20) == 0 or i == len(freq_data) - 1:
        print(f"{freq:6.0f} МГц | HPBW: {hpbw_h:5.1f}\" | "
              f"T_a_total: {res['T_a_total']:7.1f} K | "
              f"T_a_ar: {res['T_a_ar']:7.1f} K | "
              f"Вклад AR: {res['ar_contribution'] * 100:5.1f}%")

# ============================================================================
# 7. ПОДГОТОВКА ДАННЫХ ДЛЯ ГРАФИКОВ
# ============================================================================
freqs = [r['freq_MHz'] for r in results]
T_a_total = [r['T_a_total'] for r in results]
T_a_sun = [r['T_a_sun'] for r in results]
T_a_ar = [r['T_a_ar'] for r in results]

# Выберем три характерные частоты для схемы:
# 1. Низкая частота (широкий луч)
# 2. Средняя частота (средний луч)
# 3. Высокая частота (узкий луч)

low_idx = 0  # самая низкая частота
high_idx = len(freqs) - 1  # самая высокая частота
mid_idx = len(freqs) // 2  # средняя частота

characteristic_freqs = [
    (freqs[low_idx], low_idx, 'Низкая частота'),
    (freqs[mid_idx], mid_idx, 'Средняя частота'),
    (freqs[high_idx], high_idx, 'Высокая частота')
]

# ============================================================================
# 8. ПОСТРОЕНИЕ ГРАФИКОВ
# ============================================================================
print("\n" + "=" * 70)
print("ПОСТРОЕНИЕ ГРАФИКОВ")
print("=" * 70)

# Создаем фигуру с двумя графиками
fig = plt.figure(figsize=(15, 6))

# ------------------------------------------------
# ГРАФИК 1: Все компоненты антенной температуры
# ------------------------------------------------
ax1 = plt.subplot(1, 2, 1)

# Линии для антенных температур
ax1.plot(freqs, T_a_total, 'k-', linewidth=3, label='T_a_total (суммарная)')
ax1.plot(freqs, T_a_sun, color='orange', linewidth=2, label='T_a_sun (спокойное Солнце)')
ax1.plot(freqs, T_a_ar, 'r-', linewidth=2, label='T_a_ar (активная область)')

# Настройки графика
ax1.set_xlabel('Частота, МГц', fontsize=12)
ax1.set_ylabel('Антенная температура, K', fontsize=12)
ax1.set_title('Антенные температуры компонент', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=11, loc='best')
ax1.set_yscale('log')

# Добавим аннотации для характерных частот
colors = ['blue', 'green', 'purple']
for (freq, idx, label), color in zip(characteristic_freqs, colors):
    ax1.axvline(x=freq, color=color, linestyle='--', alpha=0.5, linewidth=1)
    ax1.annotate(f'{freq:.0f} МГц)', # \n({label}
                 xy=(freq, T_a_total[idx]),
                 xytext=(-60 if label == 'Высокая частота' else 10, -5 if label == 'Низкая частота' else -30),
                 textcoords='offset points',
                 fontsize=9,
                 color=color,
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# ------------------------------------------------
# ГРАФИК 2: Геометрическая схема (с заданными границами)
# ------------------------------------------------
ax2 = plt.subplot(1, 2, 2)

# Цветовая схема
SUN_COLOR = '#FFA500'  # Яркий оранжевый для Солнца
AR_COLOR = '#DC143C'  # Яркий красный для AR
BEAM_COLORS = ['#1E90FF', '#32CD32', '#8A2BE2']  # Яркие цвета для лучей

# Заданные границы графика
X_MIN, X_MAX = -1000, 1000  # По горизонтали
Y_MIN, Y_MAX = -1400, 1400  # По вертикали

# 1. Солнце (диаметр ~1920")
angles = np.linspace(0, 2 * np.pi, 300)
sun_x = R_sun_arcsec * np.cos(angles)  # Радиус 960"
sun_y = R_sun_arcsec * np.sin(angles)

# Солнце с градиентом
sun_patch = plt.Circle((0, 0), R_sun_arcsec, color=SUN_COLOR,
                       alpha=0.15, linewidth=1.5, edgecolor=SUN_COLOR,
                       label=f'Солнце (⌀≈{2 * R_sun_arcsec:.0f}")')
ax2.add_patch(sun_patch)
ax2.plot(sun_x, sun_y, color=SUN_COLOR, linewidth=1.2, alpha=0.7)

# 2. Активная область (диаметр 200") - может быть не на оси X!
ar_x = AR_POS_X + R_ar_arcsec * np.cos(angles)
ar_y = AR_POS_Y + R_ar_arcsec * np.sin(angles)

# AR с контуром и заливкой
ar_patch = plt.Circle((AR_POS_X, AR_POS_Y), R_ar_arcsec,
                      color=AR_COLOR, alpha=0.6,
                      label=f'AR (⌀={2 * R_ar_arcsec}")')
ax2.add_patch(ar_patch)
ax2.plot(ar_x, ar_y, color=AR_COLOR, linewidth=1.5, alpha=0.8)

# 3. Центры
ax2.plot(0, 0, 'o', color=SUN_COLOR, markersize=6, alpha=0.9,
         markeredgecolor='darkorange', markeredgewidth=1.2)
ax2.plot(AR_POS_X, AR_POS_Y, 'o', color=AR_COLOR, markersize=5, alpha=0.9,
         markeredgecolor='darkred', markeredgewidth=1.2)

# 4. Центр луча (всегда в (AR_POS_X, 0))
BEAM_CENTER_X = AR_POS_X
BEAM_CENTER_Y = 0.0
ax2.plot(BEAM_CENTER_X, BEAM_CENTER_Y, 'bx', markersize=4,
         markeredgewidth=0.75, label='Центр луча')

# 5. Линия от центра луча к AR (вертикальная, так как X одинаковы)
ax2.plot([BEAM_CENTER_X, BEAM_CENTER_X], [BEAM_CENTER_Y, AR_POS_Y],
         'k--', linewidth=1.0, alpha=0.5)

# 6. Эллипсы HPBW для трёх характерных частот
linewidths = [2.0, 1.6, 1.2]
alphas = [0.5, 0.6, 0.7]

for (freq, idx, label), color, lw, alpha in zip(characteristic_freqs,
                                                BEAM_COLORS,
                                                linewidths,
                                                alphas):
    res = results[idx]
    sigma_h = res['sigma_h']
    sigma_v = res['sigma_v']

    # Преобразуем sigma в HPBW
    hpbw_h = 2 * np.sqrt(2 * np.log(2)) * sigma_h
    hpbw_v = 2 * np.sqrt(2 * np.log(2)) * sigma_v

    # Эллипс HPBW с центром в точке (BEAM_CENTER_X, BEAM_CENTER_Y)
    t = np.linspace(0, 2 * np.pi, 150)
    ellipse_x = BEAM_CENTER_X + (hpbw_h / 2) * np.cos(t)
    ellipse_y = BEAM_CENTER_Y + (hpbw_v / 2) * np.sin(t)

    # Рисуем эллипс
    ax2.fill(ellipse_x, ellipse_y, color=color, alpha=alpha * 0.2)
    ax2.plot(ellipse_x, ellipse_y, color=color, linewidth=lw,
             alpha=alpha, linestyle='-')

# Настройки графика
ax2.set_xlabel('$\phi$, угл. сек', fontsize=11, fontweight='bold')
ax2.set_ylabel('θ, угл. сек', fontsize=11, fontweight='bold')
ax2.set_title(f'Геометрическая конфигурация наблюдения\nAR: ({AR_POS_X}, {AR_POS_Y}), Луч: ({BEAM_CENTER_X}, {BEAM_CENTER_Y})',
              fontsize=13, fontweight='bold', pad=12)

# УСТАНАВЛИВАЕМ ЗАДАННЫЕ ГРАНИЦЫ
ax2.set_xlim(X_MIN, X_MAX)
ax2.set_ylim(Y_MIN, Y_MAX)

# БОЛЕЕ КОНТРАСТНАЯ СЕТКА
ax2.grid(True, alpha=0.3, linewidth=0.7, color='gray', linestyle='-')

# ВАЖНО: Равный масштаб по осям!
ax2.set_aspect('equal', adjustable='box')

# ПРОРЕЖЕННЫЕ ПОДПИСИ ДЛЯ ОТМЕТОК НА ОСЯХ
ax2.set_xticks([-1000, -500, 0, 500, 1000])
ax2.set_xticklabels(['-1000', '-500', '0', '500', '1000'], fontsize=9)

ax2.set_yticks([-1400, -700, 0, 700, 1400])
ax2.set_yticklabels(['-1400', '-700', '0', '700', '1400'], fontsize=9)

# Компактная легенда
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor=SUN_COLOR, alpha=0.15, edgecolor=SUN_COLOR,
          linewidth=1, label=f'Солнце (⌀≈{2 * R_sun_arcsec:.0f}")'),
    Patch(facecolor=AR_COLOR, alpha=0.6, edgecolor=AR_COLOR,
          linewidth=1, label=f'AR (⌀={2 * R_ar_arcsec}")'),
    Line2D([0], [0], color='black', marker='x', markersize=8, lw=0,
           markeredgewidth=2, label='Центр луча'),
    Line2D([0], [0], color=BEAM_COLORS[0], lw=2.0,
           label=f'{characteristic_freqs[0][0]:.0f} МГц'),
    Line2D([0], [0], color=BEAM_COLORS[1], lw=1.6,
           label=f'{characteristic_freqs[1][0]:.0f} МГц'),
    Line2D([0], [0], color=BEAM_COLORS[2], lw=1.2,
           label=f'{characteristic_freqs[2][0]:.0f} МГц'),
]

# Легенда справа
ax2.legend(handles=legend_elements, fontsize=8, loc='center left',
           bbox_to_anchor=(0.9, 0.9), frameon=True, fancybox=True,
           framealpha=0.95, title="Обозначения:", title_fontsize=9)

# Аннотации размеров
props = dict(boxstyle='round', facecolor='white', alpha=0.95,
             edgecolor='gray', linewidth=0.6, pad=0.3)

# Диаметр Солнца (сверху)
ax2.annotate(f'⌀ Солнца ≈ {2 * R_sun_arcsec:.0f}"',
             xy=(0, R_sun_arcsec), xytext=(0, R_sun_arcsec + 150),
             ha='center', fontsize=9, fontweight='normal',
             bbox=props,
             arrowprops=dict(arrowstyle='->', color='darkorange',
                             alpha=0.6, linewidth=0.8, relpos=(0.5, 0)))

# Диаметр AR
ax2.annotate(f'⌀ AR = {2 * R_ar_arcsec}"',
             xy=(AR_POS_X, AR_POS_Y + R_ar_arcsec),
             xytext=(AR_POS_X - 600, AR_POS_Y + R_ar_arcsec + 50),
             ha='left', fontsize=8, fontweight='normal',
             bbox=props,
             arrowprops=dict(arrowstyle='->', color='darkred',
                             alpha=0.6, linewidth=0.8, relpos=(0, 0)))

# Смещение AR по X (от центра Солнца)
ax2.annotate(r'$\Delta\phi_{\mathrm{ar}} = ' + f'{AR_POS_X}"$',
             xy=(AR_POS_X/2, 0),
             xytext=(AR_POS_X/2-400, -150),
             ha='center', fontsize=8, fontweight='normal',
             bbox=props,
             arrowprops=dict(arrowstyle='->', color='black',
                             alpha=0.5, linewidth=0.6, relpos=(0.5, 0)))

# Смещение AR по Y (от центра луча) - теперь это важно!
vertical_offset = AR_POS_Y - BEAM_CENTER_Y
if abs(vertical_offset) > 10:  # Показываем только если заметное смещение
    ax2.annotate(r'$\Delta θ_{\mathrm{ar}} = ' + f'{vertical_offset:.0f}"$',
                 xy=(BEAM_CENTER_X, BEAM_CENTER_Y + vertical_offset/2),
                 xytext=(BEAM_CENTER_X - 600, BEAM_CENTER_Y + vertical_offset/2),
                 ha='left', fontsize=8, fontweight='normal',
                 bbox=props,
                 arrowprops=dict(arrowstyle='->', color='purple',
                                 alpha=0.6, linewidth=0.8, relpos=(0, 0.5)))

# Масштабная линейка (в правом нижнем углу)
scale_length = 500  # 500 угл. сек
scale_x = X_MAX - 150 - scale_length
scale_y = Y_MIN + 150

# Сама линейка
ax2.plot([scale_x, scale_x + scale_length],
         [scale_y, scale_y], 'k-', linewidth=2.5)
# Вертикальные черточки на концах
ax2.plot([scale_x, scale_x],
         [scale_y - 30, scale_y + 30], 'k-', linewidth=1.5)
ax2.plot([scale_x + scale_length, scale_x + scale_length],
         [scale_y - 30, scale_y + 30], 'k-', linewidth=1.5)
# Подпись
ax2.text(scale_x + scale_length / 2, scale_y - 70,
         f'{scale_length}"', ha='center', fontsize=8, fontweight='bold',
         bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.9))

# Добавим подписи к эллипсам (HPBW_H для каждой частоты)
y_offsets = [200, 300, 400]  # Смещаем подписи, чтобы не перекрывались
for (freq, idx, label), color, y_offset in zip(characteristic_freqs,
                                               BEAM_COLORS,
                                               y_offsets):
    res = results[idx]
    sigma_h = res['sigma_h']
    hpbw_h = 2 * np.sqrt(2 * np.log(2)) * sigma_h

    ax2.text(BEAM_CENTER_X + hpbw_h / 2 + 40, BEAM_CENTER_Y + y_offset,
             f'HPBW_H={hpbw_h:.0f}"',
             color=color, fontsize=7, ha='left', va='center',
             bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.8))

# Границы Солнца для наглядности
ax2.axvline(x=-R_sun_arcsec, color='orange', linewidth=0.6, alpha=0.3, linestyle=':', zorder=0)
ax2.axvline(x=R_sun_arcsec, color='orange', linewidth=0.6, alpha=0.3, linestyle=':', zorder=0)
ax2.axhline(y=-R_sun_arcsec, color='orange', linewidth=0.6, alpha=0.3, linestyle=':', zorder=0)
ax2.axhline(y=R_sun_arcsec, color='orange', linewidth=0.6, alpha=0.3, linestyle=':', zorder=0)
# Текст у границ Солнца
# ax2.text(R_sun_arcsec + 50, 0, 'Край\nСолнца',
#          rotation=90, va='center', fontsize=7, color='darkorange', alpha=0.6)
# ax2.text(-R_sun_arcsec - 70, 0, 'Край\nСолнца',
#          rotation=90, va='center', fontsize=7, color='darkorange', alpha=0.6)

plt.tight_layout(rect=[0, 0, 0.85, 1])  # Оставляем место для легенды

# ============================================================================
# ВЫВОД ИНФОРМАЦИИ О СООТВЕТСТВИИ ГРАНИЦ
# ============================================================================
print("\n" + "=" * 70)
print("ПРОВЕРКА ГРАНИЦ ГРАФИКА")
print("=" * 70)
print(f"Заданные границы по x: {X_MIN} ... {X_MAX}")
print(f"Заданные границы по y: {Y_MIN} ... {Y_MAX}")
print(f"Радиус Солнца: {R_sun_arcsec:.0f}\"")
print(f"Диаметр Солнца: {2 * R_sun_arcsec:.0f}\"")
print()

# Проверим, помещается ли Солнце в заданные границы
if R_sun_arcsec <= abs(X_MIN) and R_sun_arcsec <= X_MAX:
    print(f"✓ Солнце помещается по горизонтали: {R_sun_arcsec:.0f}\" <= {min(abs(X_MIN), X_MAX):.0f}\"")
else:
    print(f"⚠ Солнце не помещается по горизонтали: {R_sun_arcsec:.0f}\" > {min(abs(X_MIN), X_MAX):.0f}\"")

if R_sun_arcsec <= abs(Y_MIN) and R_sun_arcsec <= Y_MAX:
    print(f"✓ Солнце помещается по вертикали: {R_sun_arcsec:.0f}\" <= {min(abs(Y_MIN), Y_MAX):.0f}\"")
else:
    print(f"⚠ Солнце не помещается по вертикали: {R_sun_arcsec:.0f}\" > {min(abs(Y_MIN), Y_MAX):.0f}\"")

# Проверим, видны ли все эллипсы
print(f"\nВИДИМОСТЬ ЭЛЛИПСОВ HPBW:")
for freq, idx, label in characteristic_freqs:
    res = results[idx]
    sigma_h = res['sigma_h']
    sigma_v = res['sigma_v']
    hpbw_h = 2 * np.sqrt(2 * np.log(2)) * sigma_h
    hpbw_v = 2 * np.sqrt(2 * np.log(2)) * sigma_v

    right_edge = AR_POS_X + hpbw_h / 2
    left_edge = AR_POS_X - hpbw_h / 2
    top_edge = AR_POS_Y + hpbw_v / 2
    bottom_edge = AR_POS_Y - hpbw_v / 2

    visible_x = (left_edge >= X_MIN) and (right_edge <= X_MAX)
    visible_y = (bottom_edge >= Y_MIN) and (top_edge <= Y_MAX)

    if visible_x and visible_y:
        status = "✓ Полностью виден"
    elif not visible_x and visible_y:
        status = "⚠ Частично виден (выходит по горизонтали)"
    elif visible_x and not visible_y:
        status = "⚠ Частично виден (выходит по вертикали)"
    else:
        status = "✗ Не виден"

    print(f"  {freq:.0f} МГц: {status}")
sun_diameter = 2 * R_sun_arcsec

# Сохраним значение для высокой частоты
hpbw_v_high = None

for freq, idx, label in characteristic_freqs:
    res = results[idx]
    sigma_h = res['sigma_h']
    sigma_v = res['sigma_v']
    hpbw_h = 2 * np.sqrt(2 * np.log(2)) * sigma_h
    hpbw_v = 2 * np.sqrt(2 * np.log(2)) * sigma_v
    ratio = hpbw_v / sun_diameter

    # Сохраняем для высокой частоты
    if label == 'Высокая частота':
        hpbw_v_high = hpbw_v

    symbol = "✓" if hpbw_v < sun_diameter else "⚠"

    print(f"{freq:>6.0f} МГц  {hpbw_h:>10.0f}\"  {hpbw_v:>10.0f}\"  "
          f"{sun_diameter:>10.0f}\"  {ratio:>12.2f}  {symbol}")

print("-" * 75)

# Особое внимание к 3000 МГц
if hpbw_v_high is not None:
    print(f"\nДЛЯ 3000 МГЦ:")
    print(f"  Теоретический HPBW_V = 75' × (1000/3000) = 25' = 1500\"")
    print(f"  Рассчитанный HPBW_V  = {hpbw_v_high:.0f}\"")
    print(f"  Диаметр Солнца       ≈ {sun_diameter:.0f}\"")
    print(f"  HPBW_V / ⌀Солнца     = {hpbw_v_high / sun_diameter:.2f}")

    if hpbw_v_high < sun_diameter:
        print(f"  ✓ HPBW_V ({hpbw_v_high:.0f}\") МЕНЬШЕ диаметра Солнца ({sun_diameter:.0f}\")")
    else:
        print(f"  ⚠ HPBW_V ({hpbw_v_high:.0f}\") БОЛЬШЕ диаметра Солнца ({sun_diameter:.0f}\")")

    # Теоретическое значение для проверки
    hpbw_v_theoretical = 75 * 60 * (1000 / characteristic_freqs[2][0])  # 75' на 1000 МГц
    print(f"  Теоретическое значение: {hpbw_v_theoretical:.0f}\"")
    print(f"  Расхождение: {abs(hpbw_v_high - hpbw_v_theoretical):.1f}\" "
          f"({abs((hpbw_v_high - hpbw_v_theoretical) / hpbw_v_theoretical * 100):.1f}%)")
else:
    print("\nОшибка: не удалось найти данные для высокой частоты")

# Проверим, помещается ли эллипс в Солнце
print(f"\nГЕОМЕТРИЧЕСКАЯ ПРОВЕРКА:")
print(f"  Центр AR: ({AR_POS_X:.0f}, {AR_POS_Y:.0f})")
print(f"  Расстояние от центра Солнца до центра AR: {AR_POS_X:.0f}\"")

# Максимальное расстояние от центра AR до края Солнца по горизонтали
max_distance_to_edge = np.sqrt((R_sun_arcsec ** 2) - (AR_POS_Y ** 2)) - AR_POS_X
print(f"  Макс. расстояние от центра AR до края Солнца: {max_distance_to_edge:.0f}\"")

# Проверим, помещаются ли эллипсы
print(f"\nПРОВЕРКА ПОМЕЩЕНИЯ ЭЛЛИПСОВ В СОЛНЦЕ:")
for freq, idx, label in characteristic_freqs:
    res = results[idx]
    sigma_h = res['sigma_h']
    sigma_v = res['sigma_v']
    hpbw_h = 2 * np.sqrt(2 * np.log(2)) * sigma_h
    hpbw_v = 2 * np.sqrt(2 * np.log(2)) * sigma_v

    # Проверяем по горизонтали
    fits_horizontally = (AR_POS_X + hpbw_h / 2) <= (R_sun_arcsec - 50)  # С запасом 50"

    # Проверяем по вертикали
    fits_vertically = (hpbw_v / 2) <= R_sun_arcsec

    if fits_horizontally and fits_vertically:
        status = "✓ Полностью внутри Солнца"
    elif fits_vertically and not fits_horizontally:
        status = "⚠ Частично выходит за край (по горизонтали)"
    elif not fits_vertically and fits_horizontally:
        status = "⚠ Частично выходит за край (по вертикали)"
    else:
        status = "✗ Выходит за пределы Солнца"

    print(f"  {freq:.0f} МГц: HPBW_H={hpbw_h:.0f}\", HPBW_V={hpbw_v:.0f}\" - {status}")
# ============================================================================
# 10. ЭКСПОРТ РЕЗУЛЬТАТОВ
# ============================================================================
output_filename = f'sun_AR_model_D{2 * R_ar_arcsec:.0f}_pos{AR_POS_X:.0f}.csv'
with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = [
        'freq_MHz', 'HPBW_H_arcsec',
        'T_quiet_K', 'T_ar_K',
        'sigma_H_arcsec', 'sigma_V_arcsec',
        'Omega_A_arcsec2', 'I_sun_arcsec2', 'I_ar_arcsec2',
        'T_a_sun_K', 'T_a_ar_K', 'T_a_total_K',
        'ar_contribution'
    ]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for res in results:
        writer.writerow({
            'freq_MHz': res['freq_MHz'],
            'HPBW_H_arcsec': res['HPBW_H'],
            'T_quiet_K': res['T_quiet'],
            'T_ar_K': res['T_ar'],
            'sigma_H_arcsec': res['sigma_h'],
            'sigma_V_arcsec': res['sigma_v'],
            'Omega_A_arcsec2': res['Omega_A'],
            'I_sun_arcsec2': res['I_sun'],
            'I_ar_arcsec2': res['I_ar'],
            'T_a_sun_K': res['T_a_sun'],
            'T_a_ar_K': res['T_a_ar'],
            'T_a_total_K': res['T_a_total'],
            'ar_contribution': res['ar_contribution']
        })

print(f"\nРезультаты сохранены в файл: {output_filename}")

# ============================================================================
# 11. СОХРАНЕНИЕ ГРАФИКОВ В ФАЙЛ
# ============================================================================
plot_filename = f'sun_AR_model_D{2 * R_ar_arcsec:.0f}_pos{AR_POS_X:.0f}.png'
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
print(f"График сохранен в файл: {plot_filename}")

# Показать графики
plt.show()

print("\n" + "=" * 70)
print("РАСЧЁТ ЗАВЕРШЁН")
print("=" * 70)