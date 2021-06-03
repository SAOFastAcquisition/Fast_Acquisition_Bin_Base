import os
from pathlib import *


def path_to_YaDisk():

    head_path1 = 'D:\\YandexDisk'               # HP end Lenovo
    head_path2 = 'E:\\YandexDisk-svit-commerc'  # Dell
    head_path3 = 'C:\\Users\\PC\\YandexDisk'    # Stationary machine

    if os.path.isdir(head_path1):
        path = head_path1
    elif os.path.isdir(head_path2):
        path = head_path2
    elif os.path.isdir(head_path3):
        path = head_path3

    return path


if __name__ == '__main__':
    path_to_YaDisk()
    current_dir = Path.cwd()
    home_dir = Path.home()

    print(current_dir)
    print(home_dir)


