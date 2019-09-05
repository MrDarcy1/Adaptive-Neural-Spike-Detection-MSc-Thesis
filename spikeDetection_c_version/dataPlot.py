# python 3.6.8

import matplotlib.pyplot as plt

# https://blog.csdn.net/qq_35451572/article/details/79663356
import numpy as np
def fileRead(name): #read array from an file
    f = open(name)
    line = f.readline()
    data_list = []
    while line:
        num = list(map(float,line.split()))
        data_list.append(num)
        line = f.readline()
    f.close()
    return np.array(data_list)[0]


thr_file_name = r"D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/threshold_1.txt"
processed_file_name = r"D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/Processed_Data_1.txt"
data_for_thr_name = r"D:/Study/DegreeProject/Codes/spikeDetection_c_version/data/data_for_thr_1.txt"

processed_data = fileRead(processed_file_name)
threshold = fileRead(thr_file_name)
data_for_thr=fileRead(data_for_thr_name)

plt.figure(1)
plt.plot(processed_data)
plt.plot(threshold)
plt.plot(data_for_thr)

plt.show()


