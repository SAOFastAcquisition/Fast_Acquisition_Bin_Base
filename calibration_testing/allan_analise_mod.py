import numpy as np
import matplotlib.pyplot as plt
import allantools

# Генерируем данные
np.random.seed(42)
n = 5000
data = np.random.normal(0, 0.1, n) + np.cumsum(np.random.normal(0, 0.001, n))

# Различные функции анализа из allantools
sample_rate = 100  # Гц

# 1. Allan Deviation (ADEV)
taus_adev, adev, _, _ = allantools.adev(data, rate=sample_rate)

# 2. Overlapping Allan Deviation (более точная)
taus_oadev, oadev, _, _ = allantools.oadev(data, rate=sample_rate)

# 3. Modified Allan Deviation
taus_mdev, mdev, _, _ = allantools.mdev(data, rate=sample_rate)

# 4. Time Deviation
taus_tdev, tdev, _, _ = allantools.tdev(data, rate=sample_rate)

# Строим сравнительный график
plt.figure(figsize=(14, 10))

plt.subplot(2, 2, 1)
plt.loglog(taus_adev, adev, 'o-', label='ADEV')
plt.xlabel('τ')
plt.ylabel('σ(τ)')
plt.title('Allan Deviation')
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 2)
plt.loglog(taus_oadev, oadev, 's-', label='OADEV')
plt.xlabel('τ')
plt.ylabel('σ(τ)')
plt.title('Overlapping ADEV')
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 3)
plt.loglog(taus_mdev, mdev, '^-', label='MDEV')
plt.xlabel('τ')
plt.ylabel('σ(τ)')
plt.title('Modified ADEV')
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 4)
plt.loglog(taus_tdev, tdev, 'd-', label='TDEV')
plt.xlabel('τ')
plt.ylabel('σ(τ)')
plt.title('Time Deviation')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()