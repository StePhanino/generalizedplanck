import numpy as np
from scipy.interpolate import CubicSpline
from refractiveindex import RefractiveIndexMaterial
from hyperspy.signals import Signal1D
import hyperspy.api as hs
import matplotlib.pyplot as plt
import generalizedplanck as genp

def fit_signal(signal, components, gui=False, fit_lim=None, px=None, **kwargs):
    '''
    Fit a signal with multiple components.

    Parameters
    ----------
    signal : Signal1D
        Signal to be fitted.
    components : list of dictionaries
        Each dictonary is a component to be appended to the model.
        
        .. code:: python
            components2fit = [{'id_name' : 'Gaussian',
                               'name' : 'peak@1eV',
                               'active' : True,
                               'centre' : [value, bmin, bmax, free(bool)],
                               'fwhm' : [value, free(bool)]},
                              {'id_name' : 'Lorentzian',
                               'name' : 'peak@1.5eV',
                               'kwargs' : additional args to pass to 
                                          component__init__}
                              ]
            
        `id_name` and `name` are mandatory.
    fit_lim : tuple, optional
        Set a fit limit if needed in the form (left, right). The default is None.
    px : tuple, optional
        In case of a 2D navigation axis, it's the pixel indexes of the data 
        used to perform a first fitting attempt. The default is `None`.
    **kwargs :  
        Extra parameters for `fit` or `multifit`.

    Returns
    -------
    The fitted model
    
    Example
    -------
    
    >>> s = hs.data.luminescence_signal()
    >>> components = [{'id_name' : 'GaussianHF',
    >>>                'kwargs' : {},
    >>>                'name' : 'peak@3.28eV',
    >>>                'centre' : [3.28, 3.0, 3.7, True],
    >>>                'fwhm' : [0.75, False],
    >>>                'height' : [9e3],
    >>>               }]
    >>> fit_signal(s, components, fit_lim=(2.,4.5))
    '''
    
    if not isinstance(signal, Signal1D):
        raise NotImplementedError('Only implemented for Signal1D.')
        
    if not isinstance(components, list):
        components = [components]
    
    m = signal.create_model()
    bounded = False
    for comp_dict in components:
        if 'kwargs' in comp_dict.keys():
            component_kwargs = comp_dict['kwargs']
        else: 
            component_kwargs={}
        if comp_dict['id_name'] in hs.model.components1D.__all__:
                c = getattr(hs.model.components1D, comp_dict['id_name'])(**component_kwargs)
        elif comp_dict['id_name'] in genp.components.__all__:
                c = getattr(genp.components, comp_dict['id_name'])(**component_kwargs)
        else:                                                                                                                                                                                                                                                                                                                                                                                                   
            av_components = hs.model.components1D.__all__ + genp.components.__all__
            raise ValueError(f"{comp_dict['id_name']} component is not a valid"
                             "Hyperspy 1D component. Available components are: "
                             f"{','.join(av_components)}."
                             )
        if 'active' in comp_dict.keys():
            c.active = comp_dict['active']
        for param in c.parameters:
            if param.name in comp_dict.keys():
                if len(comp_dict[param.name]) == 1:
                    param.value = comp_dict[param.name][0]
                elif len(comp_dict[param.name]) == 2:
                    param.value = comp_dict[param.name][0]
                    param.free = comp_dict[param.name][-1]
                elif len(comp_dict[param.name]) == 4:
                    param.value = comp_dict[param.name][0]
                    param.free = comp_dict[param.name][-1]
                    param.bmin = comp_dict[param.name][1]
                    param.bmax = comp_dict[param.name][2]
                    bounded = True
                else:
                    text = ('Please set parameter in the following forms:'
                            '[value, bmax, bmin, free] or [value, free].')
                    raise NotImplementedError(text)
        c.name = comp_dict['name']
        m.append(c)
    m.plot()
    plt.pause(1)
    if fit_lim:
        if isinstance(fit_lim, tuple):
            left, right  = fit_lim
        elif isinstance(fit_lim, hs.roi.SpanROI):
            left = fit_lim.left
            right = fit_lim.right
        m.set_signal_range(x1=left, x2=right)
    if gui:
        m.gui()
    else:
        if signal.axes_manager.navigation_dimension == 0: 
            m.fit(bounded=bounded, **kwargs)
            m.plot(plot_components=True)
            m.print_current_values()
        else:
            #Fit on a specific spectrum to tune starting parameters
            m.axes_manager.indices = px
            m.fit(bounded=bounded, **kwargs)
            m.plot(plot_components=True)
            plt.pause(0.05) 
            choice = input('Pre-fit ok? y/n')
            if choice.lower() == 'y':
                plt.close('all')
                print('Multifit starts over the entire map')
                m.assign_current_values_to_all()
                m.multifit(bounded=bounded, iterpath='serpentine', **kwargs)
                m.plot(plot_components=True)
                m.print_current_values()
            else:
                print('Redo a pre-fit with different parameters.')
        
    return m



def get_nk(output='ev', **kwargs):
    '''
    It takes 3 possible inputs:
        
        - custom refractive index array :
            x (np.arry), x_units(string), n (complex)(np.array) 
        - constant value on the entire spectrum:
            n (float) 
        - reference to refractiveindex :
            shelf=str, book=str, page=str
    
    It returns a dictionary with x_min, x_max, x_units and the 
    interpolation function n(x).
    
    Supported x_units: um and ev.

    Parameters
    ----------
    output : str
        x axis in 'nm' on 'ev'. The default is 'ev'.
    **kwargs :
        It accepts en (in eV) or wl (in um) and complex.n (n+jk) 
        or a reference from refractive index database (shelf, book, page)

    Returns
    -------
    Dictionary with {x_extent, x_units, n(x)}.
    
    Examples
    --------
    
    >>> n1 = get_nk(n=1)
    >>> en=np.linspace(0.5,2,100)
    >>> n2_array = (np.random.uniform(0.9, 1.1, 100) + 
                    1.j * np.random.uniform(-1, 1, 100))
    >>> n2 = get_nk(x=en, x_units='ev', n=n2_array)
    >>> n3 = get_nk(shelf='main', book='GaAs', page='Papatryfonos')

    '''
    if output not in ['nm','ev']:
        raise ValueError('Output possibilities are nm or ev')
    out={}
    if  'x' in kwargs and 'x_units' in kwargs and'n' in kwargs:
        n = kwargs.get('n')
        x_units = kwargs.get('x_units')
        x = kwargs.get('x')
        if x_units in ['um','ev']:
            if x_units == output:
                x_out = x
            elif x_units == 'um':
                if output == 'nm':
                    x_out = x*1000
                elif output == 'ev':
                    x_out = np.sort(1.23984/x)
                    n = n[::-1]
            else:
                x_out = np.sort(1239.84/x)
                n = n[::-1]
            f_n = CubicSpline(x_out,n,extrapolate=False)
            extent = (x_out[0], x_out[-1])
        else:
            raise ValueError('Only um and eV supported as x axis units.')
            
    elif 'shelf' and 'book' and 'page' in kwargs:
        
        shelf = kwargs.get('shelf')
        book = kwargs.get('book')
        page = kwargs.get('page')
        refr_object = RefractiveIndexMaterial(shelf=shelf, 
                                              book=book, 
                                              page=page)
        if output == 'ev':
            x_out = np.sort(1.23984/refr_object.material.originalData['wavelength (um)'])
            n = refr_object.material.originalData['n'][::-1]
        elif output == 'nm':
            x_out = 1000*refr_object.material.originalData['wavelength (um)']
            n = refr_object.material.originalData['n']
        
        f_n = CubicSpline(x_out,n,extrapolate=False)
        extent = (x_out[0], x_out[-1])
        
    elif 'n' in kwargs and len(kwargs)==1:
        n = complex(kwargs.get('n'))
        f_n = np.vectorize(lambda x : n)
        extent = None
        
    out['x_units'] = output
    out['x_extent'] = extent
    out['function']  = f_n
    return out
    
def compare_extent(extent, ref):
    '''
    
    Parameters
    ----------
    extent : tuple
        (extent_min, extent_max).
    ref : tuple
        (ref_min, ref_max).

    Returns
    True if extent range cover ref range, False otherwise

    '''
    text = ('Extent or ref not sorted correctly.\n'
            'Each of them should be of the form : (min, max) with min < max')
    
    if extent is not None:
        if extent[1] > extent[0] and ref[1] > ref[0]:
            if  extent[0] <= ref[0] and extent[1] >= ref[1]:
                return True
            else:
                return False
        else:
            raise Exception(text)
    else:
        if ref[1] > ref[0]:
            return True
        else:
            raise Exception(text)
            
def overlap_array(array, ref, overlap='inclusion'):
    '''
    Parameters
    ----------
    array : TYPE
        DESCRIPTION.
    ref : TYPE
        DESCRIPTION.
    overlap : TYPE, optional
        DESCRIPTION. The default is 'inclusion'.

    Raises
    ------
    ValueError
        DESCRIPTION.
    Exception
        DESCRIPTION.

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    '''
    
    from hyperspy.misc.array_tools import numba_closest_index_round
    
    if not np.all(np.diff(array)>=0):
        raise ValueError('Array is not sorted!')
    elif not np.all(np.diff(ref)>=0):
        raise ValueError('Ref is not sorted!')
    
    mi1, ma1 = ref[0], ref[-1]
    mi2, ma2 = array[0], array[-1]

    if mi1 <= mi2 <= ma1:
        if ma2 < ma1:
            #Array included in Ref
            return True
        else:
            #Array is not completely included in ref
            i = numba_closest_index_round(array, np.array(ma1)) + 1
            if overlap == 'inclusion':
                return False
            elif overlap == 'partial':
                return  True, array[:i]
    elif mi1 < ma2 < ma1 and mi2 < mi1:
        #Array is not completely included in ref.
        i = numba_closest_index_round(array, np.array(mi1))
        if overlap == 'inclusion':
            return False
        elif overlap == 'partial':
            return True, array[i:]
    else:
        return False
        #raise Exception(f'Array ({mi2}->{ma2}) has no overlap with ref ({mi1}->{ma1})')


