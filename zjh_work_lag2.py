import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

savetop = '/home/laojin/my_work/lag/result4/'

data = pd.read_csv(savetop+'B_lag_para.csv')

plt.errorbar(data['logt'].values,data['B'].values,xerr = [data['logt_erl'].values,data['logt_erh'].values],
             yerr = [data['B_erl'].values,data['B_erh'].values],fmt = '.',elinewidth=1,capsize=2)
plt.xlabel('logt')
plt.ylabel('a')
plt.savefig(savetop+'A_plot_lag_para.png')
plt.close()


