#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hyperspy.api as hs
import numpy as np
import pathlib

#See example https://hyperspy.org/hyperspy-doc/current/auto_examples/create_signal/from_tabular_data_file.html


def gaas_cl():
    '''
    Create a fake CL spectrum for GaAs
    
    Returns
    -------
    hyperspy Signal1D.

    '''
    file_gaas = pathlib.Path(__file__).parents[0] / 'gaas_cl.txt'
    x, y = np.loadtxt(file_gaas, unpack=True)
    axes = [dict(axis=x, name="Energy", units="eV")]
    s = hs.signals.Signal1D(y, axes=axes)
    s.metadata.set_item("General.title", "GaAs RT CL")
    s.metadata.set_item("Signal.signal_type", 'CL')
    s.metadata.set_item("Signal.quantity", "Normalised intensity")
    return s
    

