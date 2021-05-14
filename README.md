# BN_Flow Simulator
BN_Flow simulator is written in MATLAB to simulate the biomass-nutrient dynamics with flow in 2D for BIO2020 project. We use the Marker-and-Cell (MAC) method to solve for the flow, upwind method for advection, and Cell-Centered Finite Difference (CCFD) method in space and Backward Euler in time for diffusion-reaction. We also use the operator splitting method to handle advection and diffusion-reaction separately.

See BN_Flow_Simulator.pdf for details.

## Notation
ns = scenario number
nc = case number
nd = direction of flow
- 1 = flows from bottom to top
- 2 = flows from top to bottom
- 3 = flows from left to right
- 4 = flows from right to left
ne = numer of experiments for flow-only option

## Input data in +Scenario folder
The script called S0C0.m stands for Scenario 0 Case 0. All input data are in SI units. 

## Driver
We have two drivers. Both requires inputs of ns and nc. Other input parameters are optional.
- Flow.m solves for flow and calculate the upscaled permeability. If ne = 2, we get the permeability tensor. 
```MATLAB
[K] = Flow(ns,nc,nd,ne);
```
To use this flow-only option, you need to assign the right flags in, e.g., S0C0.m. 
```MATLAB
flag.adv = 0; % turned off
flag.DR = 0;  % turned off
flag.flow = 1; % turned on
```

### Flow examples
Given some pore-scale geometry with biofilms. For flow models, we impose the parabolic inlet boundary condition and natrual outflow condition. We use the no-slip boundary conditions on the walls parallel to the flow direction. If biofilm is impermeable, i.e., G.kb = 0, then we solve the Stokes flow. If biofilm is somewhat permeable, we solve the heterogeneous Brinkman flow. We can update G.kb in, e.g., S0C0.m. 

To get the scalar permeability, run
```MATLAB
[K] = Flow(0,0,1,1);
```

To get the full permeability tensor, run
```MATLAB
[K] = Flow(0,0,1,2);
```

- BN_Flow.m solves for the coupled flow and transport model. We have options to turn on or off advection term; this can be determined in the Input data file, e.g., S0C0.m.
```MATLAB
[end_clock] = BN_Flow(ns,nc,nd); % if nd < 3, we assume that nd = 3.
```

To turn on the advection, update, e.g., S0C0.m.
```MATLAB
flag.adv = 1;
```

### Flow and transport examples
To simulate the biofilm growth by diffusion-reaction in a realistic one-pore geometry, run
```MATLAB
[end_clock] = BN_Flow(1,1,3);
```
If the last parameter nd is not assigned, then it will assume that nd = 3.

To simulate the biofilm growth by coupled flow and biomass-nutrient model, run
```MATLAB
[end_clock] = BN_Flow(2,1,3);
```


## Acknowledgements
This work was partially supported by the National Science Foundation DMS-1912938 and DMS-1522734, and by the NSF IRD plan 2019-21 for M.~Peszynska. Any opinion, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundataion.


