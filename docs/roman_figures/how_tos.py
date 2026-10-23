import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from stpsf import display_psf, display_psf_difference, roman

#### Create Roman-specific color map
roman_dark = mpl.colors.to_rgba('#00276b')
roman_med = mpl.colors.to_rgba('#70cbff')
roman_light = mpl.colors.to_rgba('#dbfeff')

roman_light_rgb = tuple(int(255 * n) for n in roman_light[:3])
roman_med_rgb = tuple(int(255 * n) for n in roman_med[:3])
roman_dark_rgb = tuple(int(255 * n) for n in roman_dark[:3])

map_dark_to_med = mpl.colors.ListedColormap(
    np.linspace(roman_dark, roman_med, 256),
    name='roman_dark_to_med')
map_dark_to_med.set_bad('r', alpha=0)

map_dark_to_light = mpl.colors.ListedColormap(
    np.linspace(roman_dark, roman_light, 256),
    name='roman_dark_to_light')
map_dark_to_light.set_bad('r', alpha=0)

#### Create stpsf_roman_page_header.png
long = 4
wide = 3

fig = plt.figure(figsize=(10, 8))
gs = fig.add_gridspec(wide, long, hspace=0.2, wspace=-0.15)
ax = gs.subplots(sharey=True, sharex=True)
axes = ax.flatten()

wfi = roman.WFI()

all_filters = [f for f in wfi.filter_list]

for i, ifilter in enumerate(sorted(all_filters)):
    ax = axes[i]

    wfi.filter = ifilter

    nlambda = None  # use defaults
    if wfi.filter in ['PRISM', 'GRISM0', 'GRISM1']:
        nlambda = 1

    psf = wfi.calc_psf(oversample=4, nlambda=nlambda)

    display_psf(psf, ax=ax, colorbar=False, title=ifilter, cmap=map_dark_to_light)

    if i not in [0, 4, 8]:
        ax.tick_params(axis='y', length=0)
    if i == 7:
        ax.tick_params(axis='x', reset=True, top=False)

    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)

axes[-1].remove()

# fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                          'stpsf_roman_page_header.png'),
#             dpi=100, facecolor='w')

#### Create compare_wfi09_wfi17.png
wfi2 = roman.WFI()
wfi2.filter = 'F129'
wfi2.detector = 'WFI09'
wfi2.detector_position = (4, 4)
psf_wfi09 = wfi2.calc_psf()
wfi2.detector = 'WFI17'
wfi2.detector_position = (4092, 4092)
psf_wfi17 = wfi2.calc_psf()

fig2, (ax_wfi09, ax_wfi17, ax_diff) = plt.subplots(1, 3, figsize=(16, 4))

display_psf(psf_wfi09, ax=ax_wfi09, imagecrop=2.0,
            title='WFI09, bottom left - F129', cmap=map_dark_to_light)
display_psf(psf_wfi17, ax=ax_wfi17, imagecrop=2.0,
            title='WFI17, top right - F129', cmap=map_dark_to_light)
display_psf_difference(psf_wfi09, psf_wfi17, ax=ax_diff, imagecrop=2.0,
                       vmax=5e-3, title='WFI09 - WFI17', cmap=map_dark_to_light)
fig2.tight_layout(w_pad=.5)
# fig2.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                          'compare_wfi09_wfi17.png'),
#              dpi=100, facecolor='w')
