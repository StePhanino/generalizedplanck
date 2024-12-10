import numpy as np
from scipy.interpolate import CubicSpline
from refractiveindex import RefractiveIndexMaterial

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


