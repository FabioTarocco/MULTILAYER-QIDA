import os
f = open("./output.ndx", "r")
all_ = set()
for line in f:
    data = line.split()
    try:
        for item in data:
            a = int(item)
            all_.add(a for a in data)
    except:
        pass
print(all_, len(all_))