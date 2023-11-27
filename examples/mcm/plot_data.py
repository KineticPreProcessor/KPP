#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2023-11-23 19:14:13 sander>

# https://pandas.pydata.org/docs/getting_started/intro_tutorials/04_plotting.html
# https://pandas.pydata.org/docs/user_guide/visualization.html

import pandas as pd
import matplotlib.pyplot as plt

# read white-space separated data files with pandas:
df_jval   = pd.read_csv('jval.dat',   index_col=0, sep='\s+')
df_mixrat = pd.read_csv('mixrat.dat', index_col=0, sep='\s+')

# print(df.head())
# print(df.columns)
# print(df.shape)

# df_mixrat['O3'] = 1E9 * df_mixrat['O3']

# change time from seconds to hours:
df_mixrat.index = df_mixrat.index / 3600
df_jval.index   = df_jval.index   / 3600

#plot = df_mixrat["O3"].plot(title="DataFrame Plot")

# plot mixing ratios:
df_mixrat.plot(figsize=(12, 7), title='mixing ratios', subplots=True)
plt.savefig('plot_data_mixrat.pdf')

# plot j-values:
df_jval.plot(figsize=(12, 7), title='j-values', subplots=True)
plt.savefig('plot_data_jval.pdf')
