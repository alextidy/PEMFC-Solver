# PEMFC-Solver
Model of the anode of a hydrogen fuel cell, in particular analysing starvation and the coupling between the upstream flowfield and electrochemical mass sink. 

Areas of improvement:
1. Current numerical implementation is explicit. Add implicit / more stable methods for mass transport to avoid extreme sensitivity to dt
2. Only models the anode. Add a corresponsing cathode model and couple them
3. Add more freedom with geometry
