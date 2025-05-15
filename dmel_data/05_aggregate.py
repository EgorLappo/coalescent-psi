#!/usr/bin/env python

import pandas as pd 

d = pd.read_parquet('~/Desktop/coalescent-psi/dmel_glinka/data/boots')

epsi_boots = d.groupby('boot').agg({'psi': 'mean'}).reset_index()['psi']
vpsi_boots = d.groupby('boot').agg({'psi': 'var'}).reset_index()['psi']

print(f'median epsi: {epsi_boots.median()} ({epsi_boots.quantile(0.025)}, {epsi_boots.quantile(0.975)})')
print(f'median vpsi: {vpsi_boots.median()} ({vpsi_boots.quantile(0.025)}, {vpsi_boots.quantile(0.975)})')