import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

m = pd.DataFrame(index=range(366), columns=range(182), dtype=float)
for delay in range(0, 182, 1):
    deaths = pd.read_csv(f'data/two_doses/deaths_{delay}.csv', header=None)
    sum_deaths = deaths.sum(axis=1)
    sum_deaths = [float(sum_deaths[0])]*(40) + list(sum_deaths)
    m.loc[delay:, delay] = sum_deaths

plt.figure(figsize=(8, 6))
ax = sns.heatmap(m, cmap='viridis_r', mask=m.isnull(), rasterized=True)
min_value = m.min().min()
min_position = np.where(m == min_value)
plt.scatter(min_position[1]+0.5, min_position[0]+0.5, color='red', s=50, label='Minimum deaths')
ax.invert_yaxis()
ax.set_yticks(np.arange(0.5, 370.5, 20))
ax.set_yticklabels(np.arange(0, 370, 20), fontsize=18, rotation=0)
ax.set_xticks(np.arange(0.5, 182.5, 20))
ax.set_xticklabels(np.arange(0, 182, 20), fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Deaths', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel(r'Start of second booster vaccination', labelpad=10, fontsize=20)
plt.xlabel(r'Start of first booster vaccination', labelpad=10, fontsize=20)
plt.tight_layout()
plt.savefig("figures/two_doses.svg", bbox_inches="tight")
plt.show()

delay_deaths = []
deaths = pd.read_csv(f'data/re1.5/deaths_2.csv', header=None, dtype=float)
delay_deaths.append(list(deaths.sum(axis=1)))
d = np.array(delay_deaths)

plt.figure(figsize=(8, 6))
plt.plot(d[0, :], label='Strategy 2', lw=2, color='purple')
plt.axhline(min_value, label='Two doses of existing vaccine', lw=2, color='red')
plt.ylabel('Deaths', fontsize=20, labelpad=10)
plt.xlabel('Start of variant-adapted vaccine deployment (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(bottom=0)
plt.xlim(left=0, right=180)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig("figures/two_doses_line.svg", bbox_inches="tight")
plt.show()
