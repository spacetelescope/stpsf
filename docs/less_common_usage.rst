*************************
Less Common Usage Options
*************************

This page documents options which exist, but are less frequently likely to be of use.  This section serves as a catch-all for some more esoteric customizations and applications. See also the :ref:`more_examples` page.


Adjusting simulated telescope line-of-sight pointing jitter
------------------------------------------------------------

Imprecisions in telescope pointing can have the effect of smearing out the PSF. WebbPSF models these as a Gaussian convolution. To simulate this with WebbPSF, the option names are ``jitter`` and ``jitter_sigma``.

>>> instrument.options['jitter'] = 'gaussian'   # jitter model name or None
>>> instrument.options['jitter_sigma'] = 0.009  # in arcsec per axis, default 0.007




Writing out only downsampled images
-----------------------------------

Perhaps you may want to calculate the PSF using oversampling, but to save disk space you only want to write out the PSF downsampled to detector resolution.

   >>> result =  inst.calc_psf(args, ...)
   >>> result['DET_SAMP'].writeto(outputfilename)

Or if you really care about writing it as a primary HDU rather than an extension, replace the 2nd line with

   >>> pyfits.PrimaryHDU(data=result['DET_SAMP'].data, header=result['DET_SAMP'].header).writeto(outputfilename)


Writing out intermediate images
-------------------------------

Your calculation may involve intermediate pupil and image planes (in fact, it most likely does). WebbPSF / POPPY allow you to inspect the intermediate pupil and image planes visually with the display keyword argument to :py:meth:`~webbpsf.JWInstrument.calc_psf`. Sometimes, however, you may want to save these arrays to FITS files for analysis. This is done with the ``save_intermediates`` keyword argument to :py:meth:`~webbpsf.JWInstrument.calc_psf`.

The intermediate wavefront planes will be written out to FITS files in the current directory, named in the format ``wavefront_plane_%03d.fits``. You can additionally specify what representation of the wavefront you want saved with the ``save_intermediates_what`` argument to :py:meth:`~webbpsf.JWInstrument.calc_psf`. This can be ``all``, ``parts``, ``amplitude``, ``phase`` or ``complex``, as defined as in :py:meth:`poppy.Wavefront.asFITS`. The default is to write ``all`` (intensity, amplitude, and phase as three 2D slices of a data cube).

If you pass ``return_intermediates=True`` as well, the return value of calc_psf is then ``psf, intermediate_wavefronts_list`` rather than the usual ``psf``.

.. warning::

   The ``save_intermediates`` keyword argument does not work when using parallelized computation, and WebbPSF will fail with an exception if you attempt to pass ``save_intermediates=True`` when running in parallel. The ``return_intermediates`` option has this same restriction.


Providing your own OPDs or pupils from some other source
--------------------------------------------------------

It is straight forward to configure an Instrument object to use a pupil OPD file of your own devising, by setting the ``pupilopd`` attribute of the Instrument object:

        >>> niriss = webbpsf.NIRISS()
        >>> niriss.pupilopd = "/path/to/your/OPD_file.fits"

If you have a pupil that is an array in memory but not saved on disk, you can pass it in as a fits.HDUList object :

        >>> myOPD = some_function_that_returns_properly_formatted_HDUList(various, function, args...)
        >>> niriss.pupilopd = myOPD

Likewise, you can set the pupil transmission file in a similar manner by setting the ``pupil`` attribute:

        >>> niriss.pupil = "/path/to/your/OPD_file.fits"


Please see the documentation for :py:class:`poppy.FITSOpticalElement` for information on the required formatting of the FITS file.
In particular, you will need to set the `PUPLSCAL` keyword, and OPD values must be given in units of meters.




Subclassing a JWInstrument to add additional functionality
----------------------------------------------------------

Perhaps you want to modify the OPD used for a given instrument, for instance to
add a defocus. You can do this by subclassing one of the existing instrument
classes to override the :py:meth:`JWInstrument._addAdditionalOptics` function. An :py:class:`OpticalSystem <poppy.OpticalSystem>` is
basically a list so it's straightforward to just add another optic there. In
this example it's a lens for defocus but you could just as easily add another
:py:class:`FITSOpticalElement <poppy.FITSOpticalElement>` instead to read in a disk file.


Note, we do this as an example here to show how to modify an instrument class by
subclassing it, which can let you add arbitrary new functionality.
There's an easier way to add defocus specifically; see below.


    >>> class FGS_with_defocus(webbpsf.FGS):
    >>>     def __init__(self, *args, **kwargs):
    >>>         webbpsf.FGS.__init__(self, *args, **kwargs)
    >>>         # modify the following as needed to get your desired defocus
    >>>         self.defocus_waves = 0
    >>>         self.defocus_lambda = 4e-6
    >>>     def _addAdditionalOptics(self, optsys, *args, **kwargs):
    >>>         optsys = webbpsf.FGS._addAdditionalOptics(self, optsys, *args, **kwargs)
    >>>         lens = poppy.ThinLens(
    >>>             name='FGS Defocus',
    >>>             nwaves=self.defocus_waves,
    >>>             reference_wavelength=self.defocus_lambda
    >>>         )
    >>>         lens.planetype = poppy.PUPIL  # tell propagation algorithm which this is
    >>>         optsys.planes.insert(1, lens)
    >>>         return optsys
    >>>
    >>> fgs2 = FGS_with_defocus()
    >>> # apply 4 waves of defocus at the wavelength
    >>> # defined by FGS_with_defocus.defocus_lambda
    >>> fgs2.defocus_waves = 4
    >>> psf = fgs2.calc_psf()
    >>> webbpsf.display_psf(psf)


Defocusing an instrument
--------------------------------

The instrument options dictionary also lets you specify an optional defocus
amount.  You can specify both the wavelength at which it should be applied, and
the number of waves of defocus (at that wavelength, specified as waves
peak-to-valley over the circumscribing circular pupil of JWST).


   >>> nircam.options['defocus_waves'] = 3.2
   >>> nircam.options['defocus_wavelength'] = 2.0e-6



