import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------------------
# Deaths

delay_deaths = []
for multiplier in [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    deaths = pd.read_csv(f'data/re3/deaths_{multiplier}.csv', header=None, dtype=float)
    delay_deaths.append(list(deaths.sum(axis=1)))

d = np.array(delay_deaths[::-1])

# Define axis labels
x_labels = ['2.5', '2.25', '2', '1.75', '1.5', '1.25', '1']

diff = -(d - min(d[-1]))
final_points = []
for i in range(len(diff)):
    for j in range(len(diff[0])):
        if diff[i][j] >= 0:
            if all(diff[i][j+k] < 0 for k in range(1, len(diff[0]) - j)):
                final_points.append(j)

plt.figure(figsize=(8, 6))
diff = diff[:-1, :]
final_points = final_points[:-1]
x_labels = x_labels[:-1]
ax = sns.heatmap(diff.T, cmap='bwr_r', center=0, vmax=10, vmin=-15)
# for each row in diff, do a vertical black line at the corresponding final_points value
for i in range(len(diff)):
    if i < len(final_points):
        ax.hlines(y=final_points[i], xmin=i, xmax=i+1, color='black', linewidth=2, linestyle='--')
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_yticks([0+0.5, 30+0.5, 60+0.5, 90+0.5, 120+0.5, 150+0.5, 180+0.5])
ax.set_yticklabels(np.linspace(0, 180, 7, dtype=int), fontsize=18)
ax.set_xticklabels(x_labels, fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Number fewer deaths compared to' '\n' 'optimal deployment of existing vaccine', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel('Start of booster vaccination (days)', labelpad=10, fontsize=20)
plt.xlabel(r'Vaccine efficacy ($\sigma$)', labelpad=10, fontsize=20)
plt.tight_layout()
# plt.savefig("figures/re3_deaths.svg", bbox_inches="tight")
plt.show()

# ---------------------------------------------------------------------------------------
# YLL

delay_yll = []
for multiplier in [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    loaded_data = pd.read_csv(f'data/re1.5/age_deaths_{multiplier}.csv', header=None, dtype=float)
    yll_ = []
    for delay in np.linspace(0, 180, 181, dtype=int):
        total_deaths = loaded_data.loc[delay, :]
        life_expectancy = pd.read_csv('data/Life_expectancy.csv', header=None, dtype=float)
        life_expectancy = life_expectancy.to_numpy().flatten().tolist()
        yll = np.multiply(total_deaths.to_numpy(), life_expectancy)
        yll_sum = np.sum(yll)
        yll_.append(yll_sum)
    delay_yll.append(yll_)

y = np.array(delay_yll[::-1])

x_labels = ['2.5', '2.25', '2', '1.75', '1.5', '1.25', '1']

final_points = []
diff = -(y - min(y[-1]))
for i in range(len(diff)):
    for j in range(len(diff[0])):
        if diff[i][j] >= 0:
            if all(diff[i][j+k] < 0 for k in range(1, len(diff[0]) - j)):
                final_points.append(j)

plt.figure(figsize=(8, 6))
diff = diff[:-1, :]
final_points = final_points[:-1]
x_labels = x_labels[:-1]
ax = sns.heatmap(diff.T, cmap='bwr_r', center=0, vmax=150, vmin=-200)
# for each row in diff, do a vertical black line at the corresponding final_points value
for i in range(len(diff)):
    if i < len(final_points):
        ax.hlines(y=final_points[i], xmin=i, xmax=i+1, color='black', linewidth=2, linestyle='--')
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_yticks([0+0.5, 30+0.5, 60+0.5, 90+0.5, 120+0.5, 150+0.5, 180+0.5])
ax.set_yticklabels(np.linspace(0, 180, 7, dtype=int), fontsize=18)
ax.set_xticklabels(x_labels, fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Number fewer YLL compared to' '\n' 'optimal deployment of existing vaccine', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel('Start of booster vaccination (days)', labelpad=10, fontsize=20)
plt.xlabel(r'Vaccine efficacy ($\sigma$)', labelpad=10, fontsize=20)
plt.tight_layout()
# plt.savefig("figures/re1.5_yll.svg", bbox_inches="tight")
plt.show()

# ---------------------------------------------------------------------------------------
# Hospitalisations

delay_hosps = []
for multiplier in [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    deaths = pd.read_csv(f'data/re3/hosp_{multiplier}.csv', header=None, dtype=float)
    delay_hosps.append(list(deaths.sum(axis=1)))

h = np.array(delay_hosps[::-1])

# Define axis labels
x_labels = ['2.5', '2.25', '2', '1.75', '1.5', '1.25', '1']


diff = -(h - min(h[-1]))
final_points = []
for i in range(len(diff)):
    for j in range(len(diff[0])):
        if diff[i][j] >= 0:
            if all(diff[i][j+k] < 0 for k in range(1, len(diff[0]) - j)):
                final_points.append(j)

plt.figure(figsize=(8, 6))
diff = diff[:-1, :]
final_points = final_points[:-1]
x_labels = x_labels[:-1]
ax = sns.heatmap(diff.T, cmap='bwr_r', center=0, vmax=60, vmin=-60)
# for each row in diff, do a vertical black line at the corresponding final_points value
for i in range(len(diff)):
    if i < len(final_points):
        ax.hlines(y=final_points[i], xmin=i, xmax=i+1, color='black', linewidth=2, linestyle='--')
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_yticks([0+0.5, 30+0.5, 60+0.5, 90+0.5, 120+0.5, 150+0.5, 180+0.5])
ax.set_yticklabels(np.linspace(0, 180, 7, dtype=int), fontsize=18)
ax.set_xticklabels(x_labels, fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Number fewer hospitalisations' '\n' 'compared to optimal deployment of' '\n' 'existing vaccine', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel('Start of booster vaccination (days)', labelpad=10, fontsize=20)
plt.xlabel(r'Vaccine efficacy ($\sigma$)', labelpad=10, fontsize=20)
plt.tight_layout()
# plt.savefig("figures/re3_hosps.svg", bbox_inches="tight")
plt.show()

# ---------------------------------------------------------------------------------------
# Cases

delay_inf = []
for multiplier in [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    deaths = pd.read_csv(f'data/re3/infections_{multiplier}.csv', header=None, dtype=float)
    delay_inf.append(list(deaths.sum(axis=1)))

i = np.array(delay_inf[::-1])

# Define axis labels
x_labels = ['2.5', '2.25', '2', '1.75', '1.5', '1.25', '1']

diff = -(i - min(i[-1]))
final_points = []
for i in range(len(diff)):
    for j in range(len(diff[0])):
        if diff[i][j] >= 0:
            if all(diff[i][j+k] < 0 for k in range(1, len(diff[0]) - j)):
                final_points.append(j)

plt.figure(figsize=(8, 6))
diff = diff[:-1, :]
final_points = final_points[:-1]
x_labels = x_labels[:-1]
ax = sns.heatmap(diff.T, cmap='bwr_r', center=0, vmax=70000, vmin=-20000)
# for each row in diff, do a vertical black line at the corresponding final_points value
for i in range(len(diff)):
    if i < len(final_points):
        ax.hlines(y=final_points[i], xmin=i, xmax=i+1, color='black', linewidth=2, linestyle='--')
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_yticks([0+0.5, 30+0.5, 60+0.5, 90+0.5, 120+0.5, 150+0.5, 180+0.5])
ax.set_yticklabels(np.linspace(0, 180, 7, dtype=int), fontsize=18)
ax.set_xticklabels(x_labels, fontsize=18)
colorbar = ax.collections[0].colorbar
colorbar.ax.tick_params(labelsize=18)
colorbar.set_label('Number fewer infections compared to' '\n' 'optimal deployment of existing vaccine', rotation=90, labelpad=10, fontsize=18)
plt.xticks(rotation=0)
plt.ylabel('Start of booster vaccination (days)', labelpad=10, fontsize=20)
plt.xlabel(r'Vaccine efficacy ($\sigma$)', labelpad=10, fontsize=20)
plt.tight_layout()
# plt.savefig("figures/re3_infs.svg", bbox_inches="tight")
plt.show()
