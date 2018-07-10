# Enthalpy formulation of convection-coupled phaseflow

The solution of the heat equation using enthalpy formulation is achieved by using FEniCS.

## Getting Started
For implementing this code, one requires a FEnics container in docker or fenics obtained through pip. A quick look on how to obtain fenics for your operating system.

### Prerequisites
FEniCS, Python Environment. 

### Installing
https://fenicsproject.org/download/

### Theory behind the code
Heat equation is given as 
                
                        𝜕T−(1/𝑆𝑡𝑒)*𝜕𝜙+∇.(𝑇𝑢)−(1/𝑃𝑟)*∆T=0
                        
Where, Ste – Stefan number = (𝐶𝑝 ∆𝑇)/𝐿   ; 
  
  𝜙 - Semi-phase field;
      
  Pr – Prandl number = (𝐶*𝑝*μ)/𝜆 ;
    
Two phase flow is interpolated into a single phase with values in the range of \[0,1\].
                                          
Where the functions are differentiated with respect to t.
            
Writing it in terms of enthalpy,

                        𝜌(𝜕ℎ + 𝑢∙∇ℎ) = 𝑘 * 𝑑𝑇/𝑑ℎ * ∆ℎ 
                   
In order to solve this equation we need to take the weak formulation, since we cannot assume the smoothness of our enthalpy function.
Multiplying by a test function v and applying Ritz Galerkin Method, we get

                        a(h,v) = −(𝑑𝑇/𝑑ℎ)*(𝑘/𝜌)*(∇𝑣,∇ℎ)
			            	      l(v) = (𝑣,𝜕𝑡ℎ)

Since this is a non-linear implementation we can write the function F as,

                        F(u,v) = ((𝑑𝑇/𝑑ℎ)(𝑘/𝜌)*(∇𝑣,∇ℎ))+(𝑣,𝜕𝑡ℎ)
                        
Solving by Newton's Method,
              
                        F(u,v) = 0
                        
#### Stefan Problem
For benchmarking we use the Stefan problem, which closely resembles the convection coupled problem. We take the flow to be stationary, 
that is the velocity is zero. Also the Prandl number is one.

## Deployment

                        python phase_with_enthalpy.py

## Plots
