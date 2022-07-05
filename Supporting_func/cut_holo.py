import numpy as np
import pandas as pd
import plotly.graph_objects as go

basedir = r'D:\YandexDisk\Amp_ML\Hologramm\Gary_Model\2022_06_25\2022_06_26\holo'
data_filename = r'\26062022.txt'
log_filename = r'\26062022.log'
RECORD_LENGTH = 200

# Поправка для сбитых часов на чувствительном комплексе (06.2022)
CORRECTION = 4.7

# Время после останоки каретки, по которому происходит усреднение
AVERAGE_TIME = pd.Timedelta(0.4, unit='s')

# .log надо чистить от следов сбоев: он должен содержать только записи длины 
# RECORD_LENGTH. Если есть записи с разной RECORD_LENGTH, надо для каждой длины
# сделать свой .log
log_data = pd.read_csv(
    basedir + log_filename,
    delim_whitespace=True,
    names=['date', 'time', 'position', 'skip'], parse_dates=[['date', 'time']]
)

# По случайности слово 'run' в строке в начале записи как раз попадает в 
# четвертую колонку, которая пустая (т.е. NaN) в остальных строках, поэтому
# используем это как признак начала очередной записи. Находим номера таких строк
# и делим данные на фреймы с отдельными записями.
record_start_indices = log_data[log_data['skip'].notna()].index.to_list()
log_records = {}
for i, ind in enumerate(record_start_indices):
    log_records[i] = log_data.date_time[ind + 1: ind + 1 + RECORD_LENGTH]
    log_records[i].index = np.arange(RECORD_LENGTH)

# Обычный формат чувствительного комплекса. Наш столбец четвертый.
data = pd.read_csv(
    basedir + data_filename,
    delim_whitespace=True,
    header=None,
    names=['timestamp', '0', '6cm', '1cm', 'skip2', 'skip3', 'skip4']
)

# Оставляем только время и данные с 1 см
data_clean = pd.concat([pd.to_datetime(data['timestamp'] 
        + 3 * 3600 + CORRECTION, unit='s'),
    data['1cm']], axis=1, keys=['date_time', '1cm'])

data_records = {}
data_records_averaged = {}
for i in range(len(log_records)):
    # Выбираем из данных куски, соотвтетствующие отдельным записям (в основном 
    # для изобразительных целей; можно обойтись и без них)
    tmin = log_records[i][0]
    tmax = log_records[i][RECORD_LENGTH - 1]
    data_records[i] = data_clean.loc[(data_clean.date_time >= tmin) 
        & (data_clean.date_time <= tmax 
        + pd.Timedelta(2, unit='s'))]

    means = np.zeros(RECORD_LENGTH)
    for j in range(0, RECORD_LENGTH):
        # Выбираем куски, соотвтетствующие отдельным шагам записи, и усредняем 
        # значения по времени
        tmin = log_records[i][j]
        tmax = tmin + AVERAGE_TIME
        selection = data_records[i].loc[(data_records[i].date_time >= tmin) 
            & ( data_records[i].date_time <= tmax)]
        means[j] = selection['1cm'].mean()

        # Для каждой записи создаем фрейм с моментами остановки каретки и 
        # соотвтствующими средними
        data_records_averaged[i] = pd.concat([
            log_records[i] + pd.Timedelta(AVERAGE_TIME / 2),
            pd.Series(means, index=np.arange(RECORD_LENGTH))],
            axis=1, keys=['date_time', '1cm'])

# Рисуем данные и средние значения - посмотреть, как они совмещаются
fig = go.Figure()
for i in range(len(data_records)):
    fig.add_trace(go.Scatter(
        x=data_records[i]['date_time'], 
        y=data_records[i]['1cm']))
    fig.add_trace(go.Scatter(
        x=data_records_averaged[i]['date_time'],
        y=data_records_averaged[i]['1cm'],
        mode='markers'))

# Пишем средние значения для каждой записи в отдельный файл
fig.show()
for i in range(len(data_records_averaged)):
    file_name = data_filename + '+' + log_filename + '.record_' + str(i).zfill(2)
    data_records_averaged[i]['1cm'].to_csv(basedir + file_name, 
        header=False,
        index=False)
