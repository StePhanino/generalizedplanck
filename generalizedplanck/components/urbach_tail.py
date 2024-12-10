from hyperspy.component import Component
import numpy as np

class UrbachTail(Component):
    r'''Urbach tail in the form
    
        .. math:: 
            \tau(E) = \tfrac{1}{2g} e^{-\left|\tfrac{E}{g}\right|}
        
        
        ============== ============== ========
         Variable        Parameter      Units  
        ============== ============== ========
         :math:`g`       Tail width     eV     
        ============== ============== ========
         
        Parameters
        ----------
        g : float
            Tail width in eV. The default is 0.015 (GaAs).
     '''

    def __init__(self, g=0.015):
        Component.__init__(self,('g'))
        self.name = 'Urbach tail component'
        self.g.value = g
        self.g.units = 'eV'
        
    def function(self,x):
        g = self.g.value
        return 1/(2*g)*np.exp(-np.abs(x/g))
