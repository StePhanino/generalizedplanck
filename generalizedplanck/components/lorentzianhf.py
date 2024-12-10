#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hyperspy Lorentzian modified version in terms of fwhm and height.
"""

import numpy as np
import dask.array as da

from hyperspy.component import _get_scaling_factor
from hyperspy._components.expression import Expression


def _estimate_lorentzian_parameters(signal, x1, x2, only_current):
    axis = signal.axes_manager.signal_axes[0]
    i1, i2 = axis.value_range_to_indices(x1, x2)
    X = axis.axis[i1:i2]

    if only_current is True:
        data = signal._get_current_data()[i1:i2]
        i = 0
        centre_shape = (1,)
    else:
        i = axis.index_in_array
        data_gi = [slice(None), ] * len(signal.data.shape)
        data_gi[axis.index_in_array] = slice(i1, i2)
        data = signal.data[tuple(data_gi)]
        centre_shape = list(data.shape)
        centre_shape[i] = 1

    cdf = np.cumsum(data,i)
    cdfnorm = cdf/np.max(cdf, i).reshape(centre_shape)

    icentre = np.argmin(abs(0.5 - cdfnorm), i)
    igamma1 = np.argmin(abs(0.75 - cdfnorm), i)
    igamma2 = np.argmin(abs(0.25 - cdfnorm), i)

    if isinstance(data, da.Array):
        icentre, igamma1, igamma2 = da.compute(icentre, igamma1, igamma2)

    centre = X[icentre]
    fwhm = (X[igamma1] - X[igamma2])
    height = data.max(i)

    return centre, height, fwhm


class LorentzianHF(Expression):

    r"""Cauchy-Lorentz distribution (a.k.a. Lorentzian function) component 
    expressed in terms of height and fwhm.

    .. math::

        f(x)=h\left[\frac{W^{2}}{4\left(x-x_{0}\right)^{2}
            +W^{2}}\right]

    ============== =============
    Variable        Parameter
    ============== =============
    :math:`h`       height
    :math:`W`       fwhm
    :math:`x_0`     centre
    ============== =============


    Parameters
    ----------
    h : float
        Height parameter. :math:`h=2A/(W\pi)` where A is the area under 
        the curve.
    fwhm : float
        Parameter corresponding to the full-width-at-half-maximum of the
        peak.
    centre : float
        Location of the peak maximum.
    **kwargs
        Extra keyword arguments are passed to the
        :class:`~.api.model.components1D.Expression` component.


    For convenience the `A` attribute can be used to get and set
    the area of the distribution.
    """

    def __init__(self, height=1., fwhm=1., centre=0., module=None, **kwargs):
        # We use `_gamma` internally to workaround the use of the `gamma`
        # function in sympy
        super().__init__(
            expression="height * (fwhm**2 / (4*(x - centre)**2 + fwhm**2))",
            name="Lorentzian",
            height=height,
            fwhm=fwhm,
            centre=centre,
            position="centre",
            module=module,
            autodoc=False,
            **kwargs)

        # Boundaries
        self.height.bmin = None
        self.height.bmax = None

        self.fwhm.bmin = None
        self.fwhm.bmax = None

        self.isbackground = False
        self.convolved = True

    def estimate_parameters(self, signal, x1, x2, only_current=False):
        """Estimate the Lorentzian by calculating the median (centre) and
        the interquartile range (fwhm).

        Note that an insufficient range will affect the accuracy of this
        method.

        Parameters
        ----------
        signal : :class:`~.api.signals.Signal1D`
        x1 : float
            Defines the left limit of the spectral range to use for the
            estimation.
        x2 : float
            Defines the right limit of the spectral range to use for the
            estimation.

        only_current : bool
            If False estimates the parameters for the full dataset.

        Returns
        -------
        bool

        Notes
        -----
        Adapted from gaussian.py and
        https://en.wikipedia.org/wiki/Cauchy_distribution

        Examples
        --------

        >>> g = hs.model.components1D.Lorentzian()
        >>> x = np.arange(-10, 10, 0.01)
        >>> data = np.zeros((32, 32, 2000))
        >>> data[:] = g.function(x).reshape((1, 1, 2000))
        >>> s = hs.signals.Signal1D(data)
        >>> s.axes_manager[-1].offset = -10
        >>> s.axes_manager[-1].scale = 0.01
        >>> g.estimate_parameters(s, -10, 10, False)
        True
        """

        super()._estimate_parameters(signal)
        axis = signal.axes_manager.signal_axes[0]
        centre, height, fwhm = _estimate_lorentzian_parameters(signal, x1, x2,
                                                              only_current)
        scaling_factor = _get_scaling_factor(signal, axis, centre)

        if only_current is True:
            self.centre.value = centre
            self.fwhm.value = fwhm
            #self.A.value = height * gamma * np.pi
            self.height.value = height
            if axis.is_binned:
                self.heigth.value /= scaling_factor
            return True
        else:
            if self.height.map is None:
                self._create_arrays()
            self.height.map['values'][:] = height
            if axis.is_binned:
                self.height.map['values'] /= scaling_factor
            self.height.map['is_set'][:] = True
            self.fwhm.map['values'][:] = fwhm
            self.fwhm.map['is_set'][:] = True
            self.centre.map['values'][:] = centre
            self.centre.map['is_set'][:] = True
            self.fetch_stored_values()
            return True

    @property
    def A(self):
        return self.fwhm.value * self.fwhm.value*np.pi/2

    @A.setter
    def A(self, value):
        self.height.value = 2*value / (self.fwhm.value * np.pi)