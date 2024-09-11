import matplotlib.pyplot as plt
from webbpsf import display_psf, roman

#### Create webbpsf-roman_page_header.png
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

# fig.savefig('webbpsf-roman_page_header.png', dpi=100, facecolor='w')

#### Create compare_wfi_sca09_sca17.png
wfi2 = roman.WFI()
wfi2.filter = 'F129'
wfi2.detector = 'SCA09'
wfi2.detector_position = (4, 4)
psf_sca09 = wfi2.calc_psf()
wfi2.detector = 'SCA17'
wfi2.detector_position = (4092, 4092)
psf_sca17 = wfi2.calc_psf()

fig2, (ax_sca09, ax_sca17, ax_diff) = plt.subplots(1, 3, figsize=(16, 4))

webbpsf.display_psf(psf_sca09, ax=ax_sca09, imagecrop=2.0,
                    title='WFI SCA09, bottom left - F129')
webbpsf.display_psf(psf_sca17, ax=ax_sca17, imagecrop=2.0,
                    title='WFI SCA17, top right - F129')
webbpsf.display_psf_difference(psf_sca09, psf_sca17, ax=ax_diff,
                               vmax=5e-3, title='SCA09 - SCA17', imagecrop=2.0)
fig2.tight_layout(w_pad=.5)
# fig2.savefig('compare_wfi_sca09_sca17.png', dpi=100, facecolor='w')


#### Create fig_coronagraph_spc_f770.png
cor = roman.RomanCoronagraph()
cor.mode = "CHARSPC_F770"

fig, ax = plt.subplots(figsize=(8,7))
mono_char_spc_psf = cor.calc_psf(nlambda=1, fov_arcsec=1.6, display=True)
# fig.savefig('fig_coronagraph_spc_f770.png', dpi=100, facecolor='w')
