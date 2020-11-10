# import scipy.io
# import numpy as np

import csv


with open('D:\\YandexDisk\Piton_Progects\\Amplifier_measure\\amp1_compression_43925MHz.csv') as File:
    reader = csv.reader(File)
    for row in reader:
        print(row)

with open('D:\\YandexDisk\Piton_Progects\\Amplifier_measure\\amp1_compression_43925MHz.csv') as File1:
    reader = csv.DictReader(File1)
    for row in reader:
        print(row['Power(dBm)'], row['S21(DB)'])