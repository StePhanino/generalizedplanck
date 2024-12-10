from hyperspy.component import Component
import numpy as np

class IdealSqrtAbsorption(Component):
    r'''Ideal absorption in the form:
        
        .. math::
            \alpha_{ideal}(E) = 
            \begin{cases}
            0 & \text{if} \quad E \le E_g \\
            a_0 \sqrt{\frac{E-E_g}{E_0-E_g}} & \text{if} \quad E>E_g
            \end{cases}
        
        
        ============== ================================= ========================
        Variable        Parameter                         Units                  
        ============== ================================= ========================
        :math:`E_g`     Band gap                          eV                     
        :math:`E_0`     Ref energy level                  eV                     
        :math:`a_0`     Abs  coefficient at :math:`E_0`   :math:`\text{cm}^{-1}` 
        ============== ================================= ========================
    
        Parameters
        ----------
        Eg : float
            Band gap energy. The default is 1.42 (GaAs).
        a0 : float
            Absorption coefficient at E0. The default is 14800 (cm-1).
            Be careful with units.
        E0 : float
            Reference energy for a0. The default is 1.6 (GaAs).
    
        Returns
        -------
        None.

    '''
    def __init__(self, Eg=1.4, a0=14800, E0=1.6):

        Component.__init__(self, ('Eg', 'a0', 'E0'))
        self.name = 'Ideal sqrt-shaped absorption coefficient'
        self.Eg.value = Eg
        self.Eg.units = 'eV'
        self.Eg.bmin = 0
        self.a0.value = a0
        self.a0.units = 'cm-1'
        self.a0.bmin = 0
        self.E0.value = E0
        self.E0.units = 'eV'
        self.E0.bmin = 0
        
    def function(self,x):
        Eg = self.Eg.value
        a0 = self.a0.value
        E0 = self.E0.value
        return np.piecewise(x, [x<=Eg, x>Eg], [lambda x : 0, 
                                               lambda x : a0*np.sqrt((x-Eg)/(E0-Eg))])
        
