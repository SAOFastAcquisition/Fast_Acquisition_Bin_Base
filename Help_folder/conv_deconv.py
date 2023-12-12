import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import signal
import matplotlib.pyplot as plt
from PIL import Image

# Загрузка изображения
filter = signal.gaussian(15, 3).reshape(1, -1)
filter = filter / np.sum(filter)
image = Image.open("image.jpg")
image = np.array(image)

# Создание фильтра Гаусса
filter = signal.gaussian(15, 3).reshape(1, -1)
b = signal.gaussian(15, 3)
filter = filter / np.sum(filter)

# Применение свертки к каналам изображения (R, G, B)
filtered_image = np.zeros_like(image)
for i in range(3):
    filtered_image[:, :, i] = signal.deconvolve2d(image[:, :, i], filter, mode="same", boundary="symm")

# Визуализация исходного и обработанного изображений
plt.subplot(1,2,1)
plt.imshow(image)
plt.title("Original Image")
plt.subplot(1,2,2)
plt.imshow(filtered_image)
plt.title("Filtered Image")
plt.show()