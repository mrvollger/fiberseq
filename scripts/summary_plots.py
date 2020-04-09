#!/usr/bin/env python

import pysam
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib
import sys
import os
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = [16, 9]
sns.set_context("poster")

PKL = sys.argv[1]
CCS = sys.argv[2]
O1 = sys.argv[3]
O2 = sys.argv[4]



#
# Accesibility plots
#
df = pd.read_pickle(PKL)
df["dist"] = df.tpl.diff(); mask = df.refName != df.refName.shift(1); df.loc[mask,"dist"] = np.nan

x = df["dist"].dropna()
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(16,16))
keep = x[x.between(-1, x.quantile(0.99999))]
sns.distplot( keep, kde=False, bins=range(int(keep.max()) +1), hist_kws={'log':True}, ax=ax1)
ax1.set_xlabel("Distance between m6A (accessible) sites")


BOT=50; TOP=200
sns.distplot( keep, kde=False, bins=range(BOT,TOP), ax=ax2); ax2.set_xlim(BOT,TOP)
heights = np.array([h.get_height() for h in ax2.patches] )
peaks, _= find_peaks(heights, distance=5, width=1)
for peak in peaks:
	ax2.text(peak+BOT, heights[peak], peak+BOT)

ax2.set_xlabel("Distance between m6A (accessible) sites")
fig.savefig(O1)


#
# Input read plots
#
ccs = pysam.AlignmentFile(CCS, check_sq=False)
nps, lengths = list(zip(*[ (rec.get_tag("np"), rec.query_length)  for rec in ccs.fetch(until_eof=True)]))
ccs_df = pd.DataFrame( {"np":nps, "length":lengths})


fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(16,16))

sns.distplot( ccs_df.np, kde=False, ax=ax1)
ax1.set_yscale('log')
ax1.set_xlabel("Number of passes per molecule")

sns.distplot( lengths, kde=False, bins=range(0,50000,50), ax=ax2); 
ax2.set_yscale('log')
ax2.set_xlabel("HiFi read length")

fig.savefig(O2)


