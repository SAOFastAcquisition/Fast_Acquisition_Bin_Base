import glob as gb
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import ast
from pathlib import Path
from Help_folder.paths_via_class import DataPaths


def sun_az_speed(_date, _az):
    """
    Функция вычисляет для конкретной даты и азимута азимутальную угловую скорость Солнца
    :param _date: дата наблюдения
    :param _az: int  азимут или список азимутов
    :return: вектор угловых скоростей Солнца в арксек/сек
    """
    _main_dir = _date[0:4]
    _data_dir = f'{_date[0:4]}_{_date[5:7]}_{_date[8:]}sun'
    _path_obj = DataPaths(_date, _data_dir, _main_dir)
    # Находим в папке первичных данных пути ко всем файлам-заголовкам
    paths = gb.glob(str(Path(_path_obj.primary_dir_path, "*.desc"))) # Set Pattern in glob() function
    if not paths:
        return 1

    # Из файлов-заголовков вынимаем словарь с параметрами наблюдения и в нем находим время кульминации 'T_OBS'
    d_mod = []
    altitude = []
    for s in paths:
        with open(s) as file:
            d = file.read()
        ind = d.find('{')
        # Словарь с параметрами наблюдения
        res_dict = ast.literal_eval(d[ind:])
        #
        time_culmination = res_dict['fits_words']['T_OBS'][-32:-6]
        try:
            dt_object1 = datetime.strptime(time_culmination, "%Y-%m-%dT%H:%M:%S.%f")
        except ValueError:
            dt_object1 = datetime.strptime(time_culmination, "%Y-%m-%dT%H:%M:%S")
            pass
        # Список времен кульминаций по азимутам как объектов datetime
        d_mod.append(dt_object1)

        # Угол вохождения Солнца при движении по дуге небосвода
        altitude0 = float(res_dict['fits_words']['ALTITUDE'])
        altitude.append(altitude0)

    #    Вычисление средней угловой скорости перемещения Солнца по небосводу между соседними азимутами
    #    в секундах угловых за секунду времени,
    #    14400 угловых секунд в 4 градусах перемещения по азимуту

    _sun_speed = []
    for i in range(len(paths) - 1):
        dt = d_mod[i + 1] - d_mod[i]
        dt = dt.total_seconds()
        # Длина дуги в градусах по линии движения по небосводу
        darc = np.sqrt((altitude[i + 1] - altitude[i]) ** 2 + 16) * 3600
        _sun_speed.append(14400 / dt)
        # _sun_speed.append(darc / dt)
    # _sun_speed.pop(3)
    # _sun_speed.pop(-4)

    # Азимуты, в которых определено время кульминации
    az = np.array([int(s[-8:-5]) for s in paths])
    # Средний азимут, в котором определена средняя скорость между соседними азимутами sun_speed
    x_az = np.array([(az[i + 1] + az[i]) / 2 for i in range(len(paths) - 1)])
    # x_az = np.delete(x_az, [3, -4])
    # Аппросимация угловой скорости Солнца полиномом 4-й степени "p"
    z = np.polyfit(x_az, np.array(_sun_speed), 6)
    p = np.poly1d(z)
    x = np.arange(-24, 24.1, 0.1)
    speed_calc = [p(z) for z in x]
    # fig, axes = plt.subplots(1, 1, figsize=(12, 12))
    # axes.plot(x_az, _sun_speed)
    # axes.plot(x, speed_calc)
    # axes.grid('on')
    # plt.show()
    return p(_az)


if __name__ == '__main__':
    date = '2024-01-04'
    main_dir = date[0:4]
    data_dir = f'{date[0:4]}_{date[5:7]}_{date[8:]}sun'
    path_obj = DataPaths(date, data_dir, main_dir)

    # sun_speed = sun_az_speed(str(Path(path_obj.primary_dir_path, "*.desc")), 24)
    sun_speed = sun_az_speed(date, 24)
    pass


