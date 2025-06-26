import matplotlib.pyplot as plt
import os
import stpsf
from stpsf import display_psf, roman

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

    display_psf(psf, ax=ax, colorbar=False, title=ifilter)

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

stpsf.display_psf(psf_wfi09, ax=ax_wfi09, imagecrop=2.0,
                    title='WFI09, bottom left - F129')
stpsf.display_psf(psf_wfi17, ax=ax_wfi17, imagecrop=2.0,
                    title='WFI17, top right - F129')
stpsf.display_psf_difference(psf_wfi09, psf_wfi17, ax=ax_diff,
                               vmax=5e-3, title='WFI09 - WFI17', imagecrop=2.0)
fig2.tight_layout(w_pad=.5)
# fig2.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                          'compare_wfi09_wfi17.png'),
#              dpi=100, facecolor='w')


#### Create fig_coronagraph_spc_f770.png
cor = roman.RomanCoronagraph()
cor.mode = "CHARSPC_F770"

fig3, ax3 = plt.subplots(figsize=(8,7))
mono_char_spc_psf = cor.calc_psf(nlambda=1, fov_arcsec=1.6, display=True)
# fig3.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                          'fig_coronagraph_spc_f770.png'),
#              dpi=100, facecolor='w')
