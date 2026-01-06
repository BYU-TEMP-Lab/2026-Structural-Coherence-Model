# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# # --- Load Data ---
# df = pd.read_csv("Not_Required_for_TC_Calcs\PythonPlots\LiF_PDF.csv", sep=',')  # your data file
# print(df.columns.tolist())

# r = df["r(A)"]
# g_LiF = df["Li-F"]

# # --- Find first few peaks (Li-F coordination shells) ---
# from scipy.signal import find_peaks
# peaks, _ = find_peaks(g_LiF, height=0.5)  # threshold can be tuned
# r_peaks = r.iloc[peaks[:3]].values  # first 3 peaks

# # --- Setup Figure ---
# fig, ax = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios':[1,1.5]})
# ax[0].set_aspect('equal')
# ax[0].set_xlim(-4, 4)
# ax[0].set_ylim(-4, 4)
# ax[0].axis('off')

# # --- Left Panel: Atomic structure around reference Li ---
# ax[0].scatter(0, 0, s=400, color='lightblue', edgecolor='k', label='Li⁺')

# colors = ['orange', 'gold', 'red']
# for i, rp in enumerate(r_peaks):
#     circle = plt.Circle((0,0), rp, color=colors[i], fill=False, lw=2, alpha=0.8)
#     ax[0].add_artist(circle)
#     ax[0].text(rp+0.1, 0.1, f'{rp:.2f} Å', fontsize=8)

# # --- Right Panel: PDF ---
# ax[1].plot(r, g_LiF, color='black', lw=1.5)
# ax[1].set_xlabel('r (Å)')
# ax[1].set_ylabel('g(r)')
# ax[1].set_title('Li–F Pair Distribution Function')

# # Highlight and connect peaks
# for rp in r_peaks:
#     gp = g_LiF.iloc[np.argmin(abs(r - rp))]
#     ax[1].plot(rp, gp, 'ro')
#     ax[1].axvline(rp, color='gray', ls=':')
#     ax[0].plot([0, rp], [0, 0], 'k:', lw=1)  # optional connecting line visually

# plt.tight_layout()
# plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# --- Load your PDF data ---
df = pd.read_csv("Not_Required_for_TC_Calcs\PythonPlots\LiF_PDF.csv", sep=None, engine="python")
print("Columns:", df.columns.tolist())  # helps confirm the exact column names

# --- Detect r and Li–F columns automatically ---
r_col = [c for c in df.columns if 'r' in c.lower()][0]
lf_col = [c for c in df.columns if 'li-f' in c.lower() or 'lif' in c.lower()][0]

r = df[r_col].values
g_LiF = df[lf_col].values

# --- Find first few Li–F peaks (coordination shells) ---
peaks, _ = find_peaks(g_LiF, height=1.1)
r_peaks = r[peaks[:3]]  # first 3 peaks

# --- Build real-space ionic structure (only on right-hand side) ---
ions = [(0, 0, 'cation')]
for i, r_shell in enumerate(r_peaks):
    n_ions = 6 + i * 4
    # Restrict to half-circle (right-hand side only)
    theta = np.linspace(-np.pi/2, np.pi/2, n_ions)
    for t in theta:
        x, y = r_shell * np.cos(t), r_shell * np.sin(t)
        ion_type = 'anion' if i % 2 == 0 else 'cation'
        ions.append((x, y, ion_type))

# --- Plot setup ---
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.2], hspace=0.05)
ax_pdf = fig.add_subplot(gs[0])
ax_real = fig.add_subplot(gs[1])

# --- Top: PDF ---
ax_pdf.plot(r, g_LiF, color='k', lw=2)
for rp in r_peaks:
    ax_pdf.axvline(rp, color='gray', ls='--', alpha=0.6)
ax_pdf.set_xlim(0, max(r))
ax_pdf.set_ylabel('g(r)')
ax_pdf.set_title('Li–F Pair Distribution Function')
ax_pdf.set_xticklabels([])

# --- Bottom: Real-space structure (half-side only) ---
for (x, y, t) in ions:
    color = 'tab:blue' if t == 'cation' else 'tab:orange'
    ax_real.add_patch(plt.Circle((x, y), 0.25, color=color, ec='k', alpha=0.9))

ax_real.set_xlim(0, max(r_peaks) + 1)
ax_real.set_ylim(-max(r_peaks) - 1, max(r_peaks) + 1)
ax_real.set_aspect('equal')
ax_real.axis('off')
ax_real.set_title('Real-space arrangement (right-hand side)')

plt.tight_layout()
plt.show()

