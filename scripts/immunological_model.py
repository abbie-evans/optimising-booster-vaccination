import numpy as np
import matplotlib.pyplot as plt

n_i0 = 1.13
n_v0 = 1.13
pi1 = -np.log(2)/35 # Rate parameter for the initial period of fast decay
pi2 = -np.log(2)/1000 # Rate parameter for the period of slow decay
ts = 75 # Timescale of switching between the fast and slow decay periods
k = 2.5 # Shape parameter of the logistic function
n50 = 0.091 # Immunity level at which 50% protection is conferred against infection by the novel variant
t = np.linspace(0, 300, 3000)

y = lambda t: n_v0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts))/(np.exp(pi1*ts)+np.exp(pi2*ts))
i = lambda t: n_i0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts))/(np.exp(pi1*ts)+np.exp(pi2*ts))
epsilon = lambda n: 1/(1+np.exp(-k*(np.log10(n/5.1)-np.log10(n50)))) # Existing vaccine
epsilon2 = lambda n: 1/(1+np.exp(-k*(np.log10(n*2/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=2
epsilon15 = lambda n: 1/(1+np.exp(-k*(np.log10(n*1.5/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=1.5
epsilon25 = lambda n: 1/(1+np.exp(-k*(np.log10(n*2.5/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=2.5
epsiloni = lambda n: 1/(1+np.exp(-k*(np.log10(n*3/5.1)-np.log10(n50)))) # Infection

z = np.zeros_like(t)

plt.figure(figsize=(8, 9))
plt.plot(t, epsilon(y(t)), color='blue', label='Vaccination/infection-induced (prior variant)', lw=2)
plt.plot(t, epsilon15(y(t)), color='orange', label=r'Vaccination-induced (novel variant, $\sigma$=1.5)', linestyle='--', lw=2)
plt.plot(t, epsilon2(y(t)), color='lightgreen', label=r'Vaccination-induced (novel variant, $\sigma$=2)', linestyle='-', lw=2)
plt.plot(t, epsilon25(y(t)), color='purple', label=r'Vaccination-induced (novel variant, $\sigma$=2.5)', linestyle='--', lw=2)
plt.plot(t, epsiloni(i(t)), color='gold', label='Infection-induced (novel variant)', lw=2)
plt.ylim(0, 1)
plt.gca().set_yticks([0, 0.25, 0.5, 0.75, 1])
plt.gca().set_yticklabels(['0', '0.25', '0.5', '0.75', '1'], fontsize=18)
plt.xticks(fontsize=18)
plt.ylabel(r'Protection against infection with' + '\n' +  r'novel variant $\left(f_j(\tau_{j,i}(t))\right)$', fontsize=20, labelpad=10)
plt.xlabel('Time following previous infection/vaccination' + '\n' + r'$\left(\tau_{j,i}(t)\right)$', fontsize=20, labelpad=10)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('figures/inf_protection.svg', bbox_inches='tight')
plt.show()

n50 = 0.021 # Immunity level at which 50% protection is conferred against hospitalisation in the novel variant outbreak

y = lambda t: n_v0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts))/(np.exp(pi1*ts)+np.exp(pi2*ts))
i = lambda t: n_i0*(np.exp(pi1*t+pi2*ts)+np.exp(pi2*t+pi1*ts)/(np.exp(pi1*ts)+np.exp(pi2*ts)))
epsilon = lambda n: 1/(1+np.exp(-k*(np.log10(n/5.1)-np.log10(n50)))) # Existing vaccine
epsilon2 = lambda n: 1/(1+np.exp(-k*(np.log10(n*2/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=2
epsilon15 = lambda n: 1/(1+np.exp(-k*(np.log10(n*1.5/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=1.5
epsilon25 = lambda n: 1/(1+np.exp(-k*(np.log10(n*2.5/5.1)-np.log10(n50)))) # Variant-adapted vaccine, sigma=2.5
epsiloni = lambda n: 1/(1+np.exp(-k*(np.log10(n*3/5.1)-np.log10(n50)))) # Infection

z = np.zeros_like(t)

plt.figure(figsize=(8, 9))
plt.plot(t, epsilon(y(t)), color='blue', label='Vaccination/infection-induced (prior variant)', lw=2)
plt.plot(t, epsilon15(y(t)), color='orange', label=r'Vaccination-induced (novel variant, $\sigma$=1.5)', linestyle='--', lw=2)
plt.plot(t, epsilon2(y(t)), color='lightgreen', label=r'Vaccination-induced (novel variant, $\sigma$=2)', linestyle='-', lw=2)
plt.plot(t, epsilon25(y(t)), color='purple', label=r'Vaccination-induced (novel variant, $\sigma$=2.5)', linestyle='--', lw=2)
plt.plot(t, epsiloni(i(t)), color='gold', label='Infection-induced (novel variant)', lw=2)
plt.ylim(0, 1)
plt.gca().set_yticks([0, 0.25, 0.5, 0.75, 1])
plt.gca().set_yticklabels(['0', '0.25', '0.5', '0.75', '1'], fontsize=18)
plt.xticks(fontsize=18)
plt.ylabel(r'Protection against hospitalisation $\left(g_j(\tau_{j,i}(t))\right)$', fontsize=20, labelpad=10)
plt.xlabel('Time following previous infection/vaccination' + '\n' + r'$\left(\tau_{j,i}(t)\right)$', fontsize=20, labelpad=10)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('figures/hosp_protection.svg', bbox_inches='tight')
plt.show()