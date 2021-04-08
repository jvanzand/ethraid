import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

trend_df = pd.read_csv('../data/trend_list.csv')
trend_df_nonan = trend_df.dropna(subset = ['dvdt', 'curv'])

print(len(trend_df), len(trend_df_nonan))

fig, ax = plt.subplots()
ax.scatter(trend_df_nonan.dvdt, trend_df_nonan.curv, s=5)
ax.set_xlim(-0.2, 0.2)
ax.set_ylim(-3.5e-5, 3.5e-5)
ax.set_xlabel(r'$\dot{\gamma}$')
ax.set_ylabel(r'$\ddot{\gamma}$')
plt.show()

