num_params=len(ansatz_params)
num_list=[]
for i in range (num_params):
    num_list.append(str(i))
num_list_non_ordinata=[]
num_list_non_ordinata=sorted(num_list, key=lambda x: x[0])
params=[0.]*num_params
for i in range (num_params):
    index=int(num_list_non_ordinata[i])
    params[index]=ansatz_params[i]
