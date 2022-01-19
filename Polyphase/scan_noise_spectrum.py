import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift
from pathlib import Path
from Supporting_func import Fig_plot as fp, align_spectrum, path_to_data


current_dir = Path.cwd()
home_dir = Path.home()
current_catalog = ' '
current_data_dir = ' '
file_path_data, head_path = path_to_data(current_catalog, current_data_dir)