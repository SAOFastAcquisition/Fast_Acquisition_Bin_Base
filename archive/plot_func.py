import matplotlib.pyplot as plt
import numpy as np
import pickle


def data_align():
    folder_align_path = r'F:\Fast_Acquisition\Alignment'
    align_file_name = r'\Align_coeff.bin'

    with open(folder_align_path + align_file_name, 'rb') as inp:
        calibration_frame_inp = pickle.load(inp)
    dat1 = calibration_frame_inp['spectrum_left1'][0]
    dat1 = dat1[-1::-1]
    n_dat1 = np.size(dat1)
    n_aver = n_dat1 / 128
    delta_freq = 7.8125 / n_aver
    freq_scale = np.array([1000 + delta_freq / 2 + i * delta_freq for i in range(n_dat1)])
    return dat1, freq_scale

def data_stocks():
    path_to_catalog = r'F:\Fast_Acquisition\2021\Results'
    path_to_file = r'\2021_03_27sun\2021-03-27_06+12_'
    path_to_data = path_to_catalog + path_to_file + 'stocks.npy'
    stocks_parameters = np.load(path_to_data, allow_pickle=True)
    # with open(path_to_data, 'rb') as inp:
    #     stocks_parameters = pickle.load(inp)
    stocks0 = stocks_parameters[0]
    stocks1 = stocks_parameters[1]
    time_ind = stocks_parameters[2]
    shape = np.shape(stocks0)
    freq = 1
    pass
    return stocks0, stocks1, time_ind


def simple_plot(inp_scale, inp_data):
    fig_main = plt.plot(inp_scale, inp_data)
    line_minor = plt.plot(inp_scale / 1000)
    # plt.plot(inp_scale, inp_data, 'r--*', inp_scale, 'g:+')  # В кавычках компактно заданы цвет, стиль линии и маркеры
    # plt.setp(fig_main, linestyle=':', color='r')
    plt.setp(fig_main, linestyle=':', color=(0, 1, 0, 0.9), marker='.', markerfacecolor='b')  # Четвертое число
    # задает прозрачность линии
    # plt.setp(line_minor, linestyle='-.')
    plt.grid()
    plt.show()
    pass


if __name__ == '__main__':
    y, x = data_align()
    simple_plot(x, y)
    data_stocks()


