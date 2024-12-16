# Generalised Planck components

A starting point for defining a series of components for hyperspy to implement the Generalised Planck law. 

In this first version we have defined 4 new components (in `/components` folder):

 - **Generalized Planck** component defined as:
   
   $$\Phi(E) = A(E) \times \frac{E^2}{4\pi^2\hbar^3c_0^2} \times 
                \frac{1}{e^{\frac{E - \Delta E_f}{k_BT}}-1} \text{,}$$
     
   where we assume, in this first version, a dielectric slab uniformly excited.
   In this case we can write the absorption in the following way:

   $$A(E) = \left[1-R(E)\right]\times\left[1-e^{-\alpha(E)d}\right] \text{,}$$
   
   where the absorption coefficient is defined by the following convolution
   
   $$\alpha_0(E) = \frac{1}{2g} \int_{-\infty}^{E-E_g} 
                \alpha_{ideal}(E-\epsilon) 
                \underbrace{e^{-\left|\tfrac{\epsilon}{g}\right|}}_\text{Urbach tail} 
                 d\epsilon  \text{.}$$
   
 - **Ideal Sqrt Absorption component** defined as:
   $\alpha_{ideal}(E) = a_0 \sqrt{\frac{E-E_g}{E_0-E_g}} \quad\text{if}\quad  E>E_g$ and 0 otherwise.

   If an analytical solution of the convolution with a tail (typically Urbach, stored in the `tail_type` attribute of the component) is available, it's stored in the `convolution_tail` method. See the [Urbach tail convolution](#Urbach-tail-convolution) section below.
 - **Urbach tail component** defined as:
   $$\tau(E) = \tfrac{1}{2g} e^{-\left|\tfrac{E}{g}\right|} \text{,}$$
   where `g` is the tail width.
   
 - **Reflectance component** - it calculates the R(E) from Fresnel coefficients.
   Different input possibilities are available for the refractive indexes inputs.


## Urbach tail convolution
In case of an ideal square-root absorption coefficient convoluted with an Urbach tail, the analytical solution is detailed in this [document](./doc/241216_analytic_conv.pdf). The [`timing_code.py`](./doc/timing_code.py) in the `doc` folder is a simple code to evaluate the time difference between analytical and numerical implementations.

## Examples

In the `examples` folder there is a jupyter-lab notebook [basic example](./examples/basic_example.ipynb)  where we fit some data with the new components.

## Installation

This first version is implemented as a stand-alone package in order to test it quickly by installing it in edit mode.

## To do/To be discussed
 - What do you think about this first implementation ?
 - In the `GeneralizedPlanck` component, the method `update_component` needs to be updated. When `GeneralizedPlanck` is appended to a model, we would have all other child components (absorption, reflectance and tail) to be appended too.
 - The `plot_components` method in the `GeneralizedPlanck` component has to be optimised with `Events` for a live update when fitting with a gui.
 - What do you think about the integration of the analytical solution of the convolution?
 - How to generalize to other absorption model? See for example [Greffet, J.-J., Bouchon, P., Brucoli, G. & Marquier, F. Light Emission by Nonequilibrium Bodies: Local Kirchhoff Law. 12 (2018)](https://doi.org/10.1103/PhysRevX.8.021008).


