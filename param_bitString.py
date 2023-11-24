
import numpy as np

N  = 10
maxCorrelationNumber = int(np.ceil(np.log2((N*(N-1))/2)))

n_parameters = 15
parameters_dict = {}
for i in range(0, n_parameters):
    parameters_dict[format(i,'b').zfill(maxCorrelationNumber)] =np.random.rand(1)

print (parameters_dict)


parameters_dict = {}
offset = 15
for i in range(0, n_parameters):
    parameters_dict[format(i + offset,'b').zfill(maxCorrelationNumber)] =np.random.rand(1)[0]

print (parameters_dict)