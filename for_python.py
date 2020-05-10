import ctypes

import numpy as np

adder = ctypes.cdll.LoadLibrary('/home/laojin/my_lat/python_c//for_python.so')

def py_cmp_func(a, b):
	print (type(a))
	print ("py_cmp_func", a[0], b[0])
	return a[0] - b[0]


a1 = np.array([1.,2.,3.,4.,5.])
b1 = np.array([2.,3.,4.,5.,6.])

a = (ctypes.c_float * 5)()
#a = (ctypes.c_float * 5)(1.,2.,3.,4.,5.)
#b = (ctypes.c_float * 5)(2.,3.,4.,5.,6.)
b = (ctypes.c_float * 5)()
c = (ctypes.c_float * 5)()#这里是生成数组

for index,val in enumerate(a1):
	a[index] = a1[index]
	b[index] = b1[index]
print(*c)
my_sum = adder.my_sum
aa = ctypes.c_float*3
print(ctypes.c_int(5))


print('aa',aa)
#my_sum.restype = ctypes.c_float#这里是设置返回值的位置。

my_sum(a,b,c,5)

print('sss',*c)
print('ggg',np.array(c))
for i in c:
	print(i)

