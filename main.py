import numpy as np
from Fig_plot import fig_plot as fp
#from Restriction import afc_correction as afc
x = np.linspace(0, 1, 100)
y = lambda x: np.sin(10 * x)
z = y(x)

# print(x, '\n', z)
# fp(z, z, x)
d = 3
A = np.vstack((x, z))

np.savetxt("foo.txt", A, header=(str(d) + 'abs'))
B = np.loadtxt('foo.txt')
f_in = open('foo.txt')
D1 = int((f_in.readline())[2])
f_in.close()

file_name0 = 'D:\\YandexDisk\\Measure\\Fast_Esquition\\24122019calibr\\20191224-1347_-14_-20-3'
pos = file_name0.find('calibr')
if not file_name0.find('sun') == -1:
    title2 = 'Sun intensity'
elif not file_name0.find('crab') == -1:
    title2 = 'Crab intensity'
elif not file_name0.find('calibr') == -1:
    title2 = 'Calibration'
    title0 = file_name0[-23:-2]
    title1 = '  ' + title0[0:4] + '.' + title0[4:6] + '.' + title0[6:8] + \
             ' chanell att=' + title0[14:17] + ' source att=' + title0[18:21]
    pass
else:
    title2 = []

title0 = file_name0[-19:-2] #+ file_name0[-1]
title1 = title0[0:4]+'.'+title0[4:6]+'.'+title0[6:8]+\
             ' time='+title0[9:11]+':'+title0[11:13]+' azimuth='+title0[14:17]

print(title1)


a = np.array([[1, 0, 1],
              [2, 3, 4],
              [0, 0, 7]])
A = 10 ** a

# d = [1, 2, 3, 4]
# d1 = np.fromiter(d)
# print(d1)

# d2 = np.asarray(d)
#
# a = np.array([ [2, 1], [2, 2], [4, 3] ])
# b = np.array([ 1, 3 ])
# total = a *b
# print(total)
# columns = (a[1:3, 0:2] != 0).sum()
# rows    = (a != 0).sum(1)
#
# print(columns, rows)

file_name0 = 'D:\\YandexDisk\\Measure\\Fast_Esquition\\26122019interference\\M28_004'
file = file_name0 + '.txt' # D:\YandexDisk\Measure\Fast_Esquition
f=open(file,"r")
lines=f.readlines()
result=[]
for x in lines[26:]:
    result.append(x.split(';')[0:2])
f.close()
result1 = [float(s) for s1 in result for s in s1]
result2 = np.reshape(result1,(len(result1)//2,2))
#Kp = afc_correction(freq)

pass
title1 = '  '+title0[0:4]+'.'+title0[4:6]+'.'+title0[6:8]+\
             ' chanell att='+title0[14:17]+' source att='+title0[19:21]