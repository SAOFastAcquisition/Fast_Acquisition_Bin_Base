import os
import glob as gb
import pandas as pd
import numpy as np
from pathlib import Path
import pickle

from Help_folder.paths_via_class import DataPaths


def load_stokes(_path='2023-10-25_05-24_stocks.npy'):
    _data = np.load(Path(_path), allow_pickle=True)
    _stokes_I = _data[0]
    _stokes_V = _data[1]
    _time_counter = _data[2]
    # plt.plot(_y0[:, 150])
    # plt.grid('both')
    # plt.show()
    return _stokes_I, _stokes_V, _time_counter


if __name__ == '__main__':
    date = '2023-12-15'
    main_dir = date[0:4]
    data_dir = f'{date[0:4]}_{date[5:7]}_{date[8:]}sun'

    path_obj = DataPaths(date, data_dir, main_dir)
    path_stokes = Path(str(path_obj.converted_data_file_path) + '_stocks.npy')
    path_stokes_base = Path(path_obj.converted_dir_path, 'stokes_base.npy')

    inVar = gb.glob(str(Path(path_obj.converted_dir_path, "*stocks.npy")))  # Set Pattern in iglob() function

    # Returning class type of variable
    # print(type(inVar))
    # # Printing list of names of all files that matched the pattern
    # print("List of the all the files in the directory having extension .py: ")
    # for py in inVar:
    #     print(py)

    for p_curr in inVar:
        data_file = os.path.basename(p_curr)
        stokes_I, stokes_V, time_count = load_stokes(p_curr)
        dict_stokes = {
            'date': data_file[0:10],
            'azimuth': data_file[13:16],
            'stokes_I': [stokes_I],
            'stokes_V': [stokes_V],
            'time_count': [time_count]
        }
        frame_stokes_curr = pd.DataFrame(dict_stokes)

        if not os.path.isfile(path_stokes_base):
            stokes_base = frame_stokes_curr
        else:
            with open(path_stokes_base, 'rb') as inp:
                stokes_base = pickle.load(inp)

        idx = stokes_base.loc[(stokes_base.date == dict_stokes['date'])
                              & (stokes_base.azimuth == dict_stokes['azimuth'])].index  #
        if not len(idx):
            stokes_base = pd.concat([stokes_base, frame_stokes_curr], axis=0, ignore_index=False)
        else:
            pass

        with open(path_stokes_base, 'wb') as out:
            pickle.dump(stokes_base, out)

    pass
