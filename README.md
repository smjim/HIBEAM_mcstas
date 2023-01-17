# HIBEAM_mcstas
Beam design and component simulation for HIBEAM nnbar experiment at ESS

Custom McStas component: "Venetian Blinds" differential neutron reflector is modeled after "Nested Reflector" concept in [X-Ray astronomy](https://imagine.gsfc.nasa.gov/educators/programs/xmm/mission/mirrors.html)
 - neutrons reflected off mirrors focused towards center of detector
 - gap in center of blade array to allow undeflected neutrons to pass through
 - possible future update: implement continuous deflection (change angle of blades throughout run) to account for distance fallen by slow neutrons due to gravity

![beam design](https://user-images.githubusercontent.com/78174712/213030309-1dfb4677-e7dc-4056-aa8c-ceba67a3f3f1.png)
current design for HIBEAM
![venetian blinds](https://user-images.githubusercontent.com/78174712/213030402-6fc815ac-5037-4668-b006-a181a084c8cd.png)
Venetian blinds model ([Desmos](https://www.desmos.com/calculator/ehkfioczjt))

HIBEAM experiment meant to improve upon 1991 (published 1994) ILL nnbar oscillation experiment
