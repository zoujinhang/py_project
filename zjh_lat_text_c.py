import ctypes
import numpy as np
from array import array

adder = ctypes.cdll.LoadLibrary('/home/laojin/my_lat/python_c/spectrum_tool.so')
a = np.arange(10000)*1.0
spe = (ctypes.c_double * 10000)(*list(a))
n_ret = 10000/5
print(n_ret)
ret = (ctypes.c_double* int(n_ret))()

#for index,val in enumerate(a):
#	spe[index] = val
#ret = ctypes.POINTER(ctypes.c_float)
#adder.A_spec.restype = ctypes.POINTER(ctypes.c_float)
arr = array('l', range(10))

a2 = ctypes.POINTER(ctypes.c_float)
adder.A_spec(spe,ret,10000,5,int(n_ret))

#ctypes.memmove(ret,a,ctypes.sizeof(ctypes.c_float)*int(n_ret))
print(arr)
print(len(arr) * arr.itemsize )
print(ctypes.sizeof(ctypes.c_float)*int(n_ret))
print(np.array(ret))
print(a2)
#get_value=ctypes.cast(a2, ctypes.c_float).value
#print(get_value)









