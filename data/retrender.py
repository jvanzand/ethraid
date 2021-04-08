import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

hires_resid = pd.read_csv('resid_HIRES.txt')
apf_resid = pd.read_csv('resid_APF.txt')

trend = 0.097

hires_resid['tel'] = ['hires_j' for i in range(len(hires_resid))]
apf_resid['tel'] = ['apf' for i in range(len(apf_resid))]

hires_resid.columns = ['bjd', 'mnvel', 'errvel', 'tel']
apf_resid.columns = ['bjd', 'mnvel', 'errvel', 'tel']

hires_resid['mnvel'] = hires_resid['mnvel'] + trend*(hires_resid['bjd'] - hires_resid['bjd'][0])
apf_resid['mnvel'] = apf_resid['mnvel'] + trend*(apf_resid['bjd'] - apf_resid['bjd'][0])

plt.scatter(hires_resid['bjd'], hires_resid['mnvel'], s=2)
plt.scatter(apf_resid['bjd'], apf_resid['mnvel'], s=2)

full_residuals = hires_resid.append(apf_resid).sort_values(by='bjd').reset_index(drop=True)

full_residuals.to_csv('residuals.csv', index=False)