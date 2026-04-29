import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

M = pd.read_csv('data/UK_Contacts.csv', header=None)

plt.figure(figsize=(8, 6))
sns.heatmap(M.T, cmap="Blues", vmax=4)
colorbar = plt.gca().collections[0].colorbar
colorbar.ax.tick_params(labelsize=16)
colorbar.set_ticks([0, 1, 2, 3, 4])
colorbar.set_ticklabels(['0', '1', '2', '3', '4+'])
colorbar.set_label('Number of contacts per day', rotation=90, labelpad=10, fontsize=18)
plt.gca().invert_yaxis()
plt.xlabel(r'Age group of individual ($b$)', fontsize=20, labelpad=10)
plt.ylabel(r'Age group of contact ($a$)', fontsize=20, labelpad=10)
plt.xticks(np.arange(0, 18, 2), fontsize=18)
plt.yticks(np.arange(0, 18, 2), fontsize=18)
plt.gca().set_xticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80])
plt.gca().set_yticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80])
plt.tight_layout()
plt.savefig("figures/contact_matrix.svg", bbox_inches="tight")
plt.show()