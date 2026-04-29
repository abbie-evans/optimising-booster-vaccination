import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------------------
# Deaths

delay_deaths = []
for multiplier in [1, 2]:
    deaths = pd.read_csv(f'data/re1.5/deaths_{multiplier}.csv', header=None, dtype=float)
    delay_deaths.append(list(deaths.sum(axis=1)))

d = np.array(delay_deaths[::-1])

plt.figure(figsize=(8, 6))
plt.plot([np.min(d[-1])]*len(d[0, :]), label='Strategy 1', lw=2, color='gold') # Optimal timing of existing vaccine
plt.plot(d[0, :], label='Strategy 2', lw=2, color='purple')

delay_deaths = []
optim_day = []
sens_deaths = []
sens_deaths_0_5 = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/s3/re1.5/deaths_{delay}.csv', header=None, dtype=float)
    deaths_0_5 = pd.read_csv(f'data/s3/re1.5/0.5/deaths_{delay}.csv', header=None, dtype=float)
    deaths_1_5 = pd.read_csv(f'data/s3/re1.5/1.5/deaths_{delay}.csv', header=None, dtype=float)
    min_deaths = deaths.sum(axis=1).min()
    min_deaths_1_5 = deaths_1_5.sum(axis=1).min()
    min_deaths_0_5 = deaths_0_5.sum(axis=1).min()
    optim_day.append(deaths.sum(axis=1).idxmin())
    delay_deaths.append(min_deaths)
    sens_deaths.append(min_deaths_1_5)
    sens_deaths_0_5.append(min_deaths_0_5)

plt.plot(delay_deaths, label='Strategy 3', lw=2, color='green')
# Can also plot the sensitivity analyses where different amounts of the existing vaccine is used
# plt.plot(sens_deaths_0_5, label='Strategy 3 (x0.5)', lw=2, color='green')
# plt.plot(sens_deaths, label='Strategy 3 (x1.5)', lw=2, color='green')

delay_deaths = []
sens_deaths = []
sens_deaths_0_5 = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/s4/re1.5/deaths_{delay}.csv', header=None, dtype=float)
    deaths_0_5 = pd.read_csv(f'data/s4/re1.5/0.5/deaths_{delay}.csv', header=None, dtype=float)
    deaths_1_5 = pd.read_csv(f'data/s4/re1.5/1.5/deaths_{delay}.csv', header=None, dtype=float)
    min_deaths = deaths.sum(axis=1).min()
    min_deaths_0_5 = deaths_0_5.sum(axis=1).min()
    min_deaths_1_5 = deaths_1_5.sum(axis=1).min()
    delay_deaths.append(min_deaths)
    sens_deaths.append(min_deaths_1_5)
    sens_deaths_0_5.append(min_deaths_0_5)

plt.plot(delay_deaths, label='Strategy 4', lw=2, color='dodgerblue')
# Can also plot the sensitivity analyses where different amounts of the existing vaccine is used
# plt.plot(sens_deaths_0_5, label='Strategy 4 (x0.5)', lw=2, color='dodgerblue')
# plt.plot(sens_deaths, label='Strategy 4 (x1.5)', lw=2, color='dodgerblue')

plt.ylabel('Deaths (per 100K)', fontsize=20, labelpad=10)
plt.xlabel('Start of variant-adapted vaccine deployment (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(bottom=0)
plt.xlim(left=0, right=180)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.legend(fontsize=16)
plt.savefig('figures/deaths_lines.svg', bbox_inches='tight')
plt.show()

# ---------------------------------------------------------------------------------------
# YLLs

plt.figure(figsize=(8, 6))

delay_yll = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/re1.5/age_deaths_1.csv', header=None, dtype=float).iloc[delay, :]
    life_expectancy = pd.read_csv('data/Life_expectancy.csv', header=None, dtype=float)
    life_expectancy = life_expectancy.to_numpy().flatten().tolist()
    yll = np.multiply(deaths, life_expectancy)
    delay_yll.append(np.sum(yll))

min_yll = np.min(delay_yll)
plt.plot([min_yll]*len(delay_yll), label='Strategy 1', lw=2, color='gold')

delay_yll = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/re1.5/age_deaths_2.csv', header=None, dtype=float).iloc[delay, :]
    life_expectancy = pd.read_csv('data/Life_expectancy.csv', header=None, dtype=float)
    life_expectancy = life_expectancy.to_numpy().flatten().tolist()
    yll = np.multiply(deaths, life_expectancy)
    delay_yll.append(np.sum(yll))

plt.plot(np.array(delay_yll).flatten(), label='Strategy 2', lw=2, color='purple')

delay_yll = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/s3/re1.5/age_deaths_{delay}.csv', header=None, dtype=float)
    life_expectancy = pd.read_csv('data/Life_expectancy.csv', header=None, dtype=float)
    life_expectancy = life_expectancy.to_numpy().flatten().tolist()
    yll = np.multiply(deaths, life_expectancy)
    min_yll = yll.sum(axis=1).min()
    delay_yll.append(min_yll)

plt.plot(np.array(delay_yll).flatten(), label='Strategy 3', lw=2, color='green')

delay_yll = []
for delay in np.linspace(0, 180, 181, dtype=int):
    deaths = pd.read_csv(f'data/s4/re1.5/age_deaths_{delay}.csv', header=None, dtype=float)
    life_expectancy = pd.read_csv('data/Life_expectancy.csv', header=None, dtype=float)
    life_expectancy = life_expectancy.to_numpy().flatten().tolist()
    yll = np.multiply(deaths, life_expectancy)
    min_yll = yll.sum(axis=1).min()
    delay_yll.append(min_yll)

plt.plot(np.array(delay_yll).flatten(), label='Strategy 4', lw=2, color='dodgerblue')

plt.ylabel('YLL (per 100K)', fontsize=20, labelpad=10)
plt.xlabel('Start of variant-adapted vaccine deployment (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(bottom=0)
plt.xlim(left=0, right=180)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.legend(fontsize=16)
plt.savefig('figures/yll-lines.svg', bbox_inches='tight')
plt.show()
