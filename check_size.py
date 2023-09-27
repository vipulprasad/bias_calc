import pandas as pd
import os

path = "/home/vipul/vipul/halo_clutering/bias_calc/box7500/z3.0/halo_cat/"
halo_cats = os.listdir(path)

for file in halo_cats:
    print(file)
    data = pd.read_csv(path+file)
    print(len(data))
