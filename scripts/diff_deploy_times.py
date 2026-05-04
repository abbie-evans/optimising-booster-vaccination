import matplotlib.pyplot as plt
import pandas as pd

inf_no_vacc = pd.read_csv('data/infections_400.csv', header=None, dtype=float).values.tolist()
inf_0 = pd.read_csv('data/infections_0.csv', header=None, dtype=float).values.tolist()
inf_109 = pd.read_csv('data/infections_109.csv', header=None, dtype=float).values.tolist()

death_no_vacc = pd.read_csv('data/deaths_400.csv', header=None, dtype=float).values.tolist()
death_0 = pd.read_csv('data/deaths_0.csv', header=None, dtype=float).values.tolist()
death_109 = pd.read_csv('data/deaths_109.csv', header=None, dtype=float).values.tolist()

plt.figure(figsize=(8, 6))
plt.plot(inf_no_vacc[0], label='Unmitigated', lw=2)
plt.plot(inf_0[0], label='No delay', lw=2)
plt.plot(inf_109[0], label='109 day delay', lw=2)
plt.ylabel('New infections', fontsize=20, labelpad=10)
plt.xlabel('Time (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('figures/new_infections.svg', bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
cum_inf = [sum(inf_no_vacc[0][:i+1]) for i in range(len(inf_no_vacc[0]))]
plt.plot(cum_inf, label='Unmitigated', lw=2)
cum_inf = [sum(inf_0[0][:i+1]) for i in range(len(inf_0[0]))]
plt.plot(cum_inf, label='No delay', lw=2)
cum_inf = [sum(inf_109[0][:i+1]) for i in range(len(inf_109[0]))]
plt.plot(cum_inf, label='109 day delay', lw=2)
plt.ylabel('Total infections', fontsize=20, labelpad=10)
plt.xlabel('Time (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('figures/total_infections.svg', bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(death_no_vacc[0], label='Unmitigated', lw=2)
plt.plot(death_0[0], label='No delay', lw=2)
plt.plot(death_109[0], label='109 day delay', lw=2)
plt.ylabel('New deaths', fontsize=20, labelpad=10)
plt.xlabel('Time (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('figures/new_deaths.svg', bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
cum_inf = [sum(death_no_vacc[0][:i+1]) for i in range(len(death_no_vacc[0]))]
plt.plot(cum_inf, label='Unmitigated', lw=2)
cum_inf = [sum(death_0[0][:i+1]) for i in range(len(death_0[0]))]
plt.plot(cum_inf, label='No delay', lw=2)
cum_inf = [sum(death_109[0][:i+1]) for i in range(len(death_109[0]))]
plt.plot(cum_inf, label='109 day delay', lw=2)
plt.ylabel('Total deaths', fontsize=20, labelpad=10)
plt.xlabel('Time (days)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('figures/total_deaths.svg', bbox_inches='tight')
plt.show()