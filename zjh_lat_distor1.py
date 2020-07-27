
import numpy as np
import matplotlib.pyplot as plt
import os

savedir = '/home/laojin/my_lat/distor/'
if os.path.exists(savedir) == False:
	os.makedirs(savedir)

rate_max = 5000
time_range = [0,500]

float_percentage = 0.3
float_period = 50
float_model = np.sin

w = 2*np.pi / float_period

model_t = np.linspace(time_range[0],time_range[-1],1000)

model_rate = rate_max*float_percentage*0.5*(float_model(w*model_t)+1)+rate_max-rate_max*float_percentage

t_len = time_range[-1] - time_range[0]

data1_x = np.random.rand(t_len*rate_max) * t_len - time_range[0]
data1_x = np.sort(data1_x)
data1_y = np.random.rand(t_len*rate_max) #[0,1]

data1_p = float_model(w*data1_x)
data1_p = 1 - float_percentage + (data1_p - data1_p.min())/(data1_p.max()-data1_p.min()) * float_percentage

data1_y_p = data1_p - data1_y

index_ = np.where(data1_y_p>=0)[0]

data_t = data1_x[index_]
data_y = np.random.rand(data_t.size)*500

bins_1 = np.arange(time_range[0],time_range[-1],1)
bin1_n = np.histogram(data_t,bins = bins_1)[0]
bin1_rate = bin1_n/1
bin1_rate = np.concatenate((bin1_rate[:1],bin1_rate))

bins_01 = np.arange(time_range[0],time_range[-1],0.1)
bin_01_n = np.histogram(data_t,bins = bins_01)[0]
bin_01_rate = bin_01_n/0.1
bin_01_rate = np.concatenate((bin_01_rate[:1],bin_01_rate))

bins_001 = np.arange(time_range[0],time_range[-1],0.01)
bin_001_n = np.histogram(data_t,bins = bins_001)[0]
bin_001_rate = bin_001_n/0.01
bin_001_rate = np.concatenate((bin_001_rate[:1],bin_001_rate))

bins_0001 = np.arange(time_range[0],time_range[-1],0.001)
bin_0001_n = np.histogram(data_t,bins = bins_0001)[0]
bin_0001_rate = bin_0001_n/0.001
bin_0001_rate = np.concatenate((bin_0001_rate[:1],bin_0001_rate))

plt.figure(figsize = (15,10))
plt.subplot(2,1,1)

plt.plot(data_t,data_y,',',color = 'k')
plt.xlim(time_range)
plt.yticks([])
plt.ylabel('samples')
plt.ylim(0,100)

plt.subplot(2,1,2)
plt.step(bins_0001,bin_0001_rate,label = 'bins = 0.001s')
plt.step(bins_001,bin_001_rate,label = 'bins = 0.01s')
plt.step(bins_01,bin_01_rate,label = 'bins = 0.1s')
plt.step(bins_1,bin1_rate,label = 'bins = 1s')
plt.plot(model_t,model_rate,color = 'r',label = 'model rate')
plt.xlim(time_range)
plt.ylim(0)
plt.xlabel('time')
plt.ylabel('rate')
plt.legend()
plt.savefig(savedir + 'A_samples.png')
plt.close()











