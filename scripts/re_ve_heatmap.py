import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

m = pd.DataFrame()

for i, ve in enumerate([1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]):
    for j, re in enumerate(np.arange(1.5, 3.1, 0.1)):
        deaths = pd.read_csv(f'data/optim_timing/{ve}/deaths_{np.round(re, 1)}.csv', header=None, dtype=float)
        sum_deaths = deaths.sum(axis=1)
        optim_time = sum_deaths.idxmin()
        m.loc[j, i] = optim_time

plt.figure(figsize=(8, 6))
x_labels = ['1', '1.25', '1.5', '1.75', '2', '2.25', '2.5']
ax = sns.heatmap(m, cmap='coolwarm')
ax.invert_yaxis()
ax.set_yticks(np.arange(0+0.5, 16+0.5, 1))
ax.set_yticklabels(np.round(np.linspace(1.5, 3, 16), 1), fontsize=18, rotation=0)
ax.set_xticklabels(x_labels, fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Optimal timing of vaccination', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel(r'Effective reproduction number ($R_e$)', labelpad=10, fontsize=20)
plt.xlabel(r'Vaccine efficacy ($\sigma$)', labelpad=10, fontsize=20)
plt.tight_layout()
plt.savefig("figures/re_ve.svg", bbox_inches="tight")
plt.show()
