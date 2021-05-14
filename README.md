# BN_Flow Simulator
BN_Flow simulator is written in MATLAB to simulate the biomass-nutrient dynamics with flow in 2D for BIO2020 project. We use the Marker-and-Cell (MAC) method to solve for the flow, upwind method for advection, and Cell-Centered Finite Difference (CCFD) method in space and Backward Euler in time for diffusion-reaction. We also use the operator splitting method to handle advection and diffusion-reaction separately.

See BN_Flow_Simulator.pdf for details.

## Example: Solve for the flow
```MATLAB
if 123; else 344; end
```

## Acknowledgements
This work was partially supported by the National Science Foundation DMS-1912938 and DMS-1522734, and by the NSF IRD plan 2019-21 for M.~Peszynska. Any opinion, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundataion.


