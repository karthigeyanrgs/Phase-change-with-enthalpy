# Enthalpy formulation of convection-coupled phaseflow

The solution of the heat equation using enthalpy formulation is achieved by using FEniCS.

## Getting Started
The primary motive is to be able to monitor the phase change interfaces for a phase changing liquid. Phase changing liquid comprises of a broad class of uprising materials which can be extensively used for thermal energy storage. A driving force on these phase interfaces is the convective heat transfer in the fluid. We consider the model of ice encapsulated by a heat plate acting vertically onto one segment of the ice. In order to solve this problem numerically, we employ the single domain semi-phase field approach.

### Prerequisites
FEniCS, Python Environment. 

### Installing
https://fenicsproject.org/download/

### Theoratical background.
Invoking the incompressible Navier stokes equation, but diverting our focus to the energy equation.
Energy equation is given as 
                
                        𝜕T−(1/𝑆𝑡𝑒)*𝜕𝜙+∇.(𝑇𝑢)−(1/𝑃𝑟)*∆T=0
                        
Where, Ste – Stefan number = (𝐶𝑝 ∆𝑇)/𝐿   ; 
  
  𝜙 - Semi-phase field;
      
  Pr – Prandl number = (𝐶*𝑝*μ)/𝜆 ;
     
  T - Temperature field of the substance (Here ice.)
    
  L - Latent heat of fusion. 
  
  Cp - Specific heat capacity with constant pressure.
    
Two phase flow is interpolated into a single phase with values in the range of \[0,1\].
                                          
Where the functions are differentiated with respect to t.
            
Writing it in terms of enthalpy,

                        𝜌(𝜕ℎ + 𝑢∙∇ℎ) = 𝑘 * 𝑑𝑇/𝑑ℎ * ∆ℎ
			
Here ,  𝜌 - Density of the substance (ice)
	u - Velocity of the fluid(ice - water mixture)
        k - Thermal Conductivity.
	h - Enthalpy field.
In order to solve this equation we need to take the weak formulation, since we cannot assume the smoothness of our enthalpy function.
Multiplying by a test function v and applying Ritz Galerkin Method, we get

                        a(h,v) = −(𝑑𝑇/𝑑ℎ)*(𝑘/𝜌)*(∇𝑣,∇ℎ)
			            	      l(v) = (𝑣,𝜕𝑡ℎ)
		  
Here a(h,v) is a bilinear function maps the linear polynomials of test function to the functional space which here is just an one dimensional grid , while the l(v) is a linear function in the functional space.

Since this is a non-linear implementation we can write the function F as,

                        F(u,v) = ((𝑑𝑇/𝑑ℎ)(𝑘/𝜌)*(∇𝑣,∇ℎ))+(𝑣,𝜕𝑡ℎ)
                        
Solving by Newton's Method,
              
                        F(u,v) = 0
                        
#### Stefan Problem
For benchmarking we use the Stefan problem, which closely resembles the convection coupled problem. We take the flow to be stationary, 
that is the velocity is zero. Also the Prandl number is one.

## Deployment

                        python phase_with_enthalpy.py

## Author
* **Karthigeyan Ganesh Shankar**  - [Karthigeyan](https://github.com/karthigeyanrgs)

## References
[The fenics project](https://doi.org/10.11588/ans.2015.100.20553)
[Monolithic simulation of convection-coupled phase-change -verification and reproducibility](https://arxiv.org/abs/1801.03429)
[Geo-fluid-dynamics]()
##Acknowledgement
I would like to thank Professor Julia Kowalski and Mr. Alexander G. Zimmerman 
