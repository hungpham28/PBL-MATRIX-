"""
solve 
create a system of linear equations (n variable) and create 1 new file 
    where n : level matrix
          max: limited Bigdata    
"""
from random import random
path=""
path_w = path+input("Nhập tên file :")
import random
n=int(input("Nhập n cho ma trận: "))
max=int(input("Nhập giới hạn dữ liệu cho ma trận: "))
with open(path_w, mode='w') as f:
    f.write(str(n)+"\n")
    for _ in range(n):
        s=""
        for i in range(n+1):
            s+=str(random.randint(-max,max))+" "
        f.write(s+"\n")
