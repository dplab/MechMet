# MechMet Introduction
MatLab code for computing mechanical metrics of polycrystals from 3D volume meshed data, using the same datatypes as FEpX (https://fepx.info). MechMet computes elastic deformations of virtual polycrystals for several possible loading modes such a uniaxial or biaxial extension.  With the computed responses, MechMet evaluates stiffness-based metrics based on the single-crystal elastic moduli and mechanical interactions among the crystals. Further, using the single-crystal yield surface (determined from input of slip system strengths), MechMet evaluates metrics derived from combinations of strength and stiffness that can provide insight to where yielding may occur or where damage may accumulate. Neper (https://neper.info/) can be used to create synthetic or experimentally derived 3D datasets via tessellation as initial input to MechMet.

The MechMet code was built in the Cornell Deformation Lab primarily by Paul Dawson. A forthcoming paper has been submitted, which details the use of the MechMet code, some examples, and the metrics described here. The manuscript DOI and citable text will be updated and listed here upon acceptance of the manuscript.

_"Mechanical Metrics of Virtual Polycrystals (MechMet), P.R. Dawson, M.P. Echlin, M.P. Miller, T.M. Pollock, J. Wendorf, L.H. Mills, J.C. Stinville, M.A. Charpagne, submitted, 2021"_

MechMet computes the following metrics:

#### Stiffness-based metrics:

* Directional Stiffness = 1 / S<sub>11</sub> (e.g. in the x<sub>1</sub> direction)  
* Embedded Stiffness = E&#770; = &epsilon;<sub>11</sub> / &sigma;&#773;<sub>11</sub> , where &sigma;&#773;<sub>11</sub> is the average value of &sigma;<sub>11</sub> over the polycrystal  
* Single-Crystal Directional Stiffness = (1 / S<sub>11</sub>)    
* Embedded Stiffness to the Single-Crystal Directional Stiffness (RE2SX) = E&#770; / (1 / S<sub>11</sub>)

#### Strength+Stiffness-based metrics
* Relative Schmid Factor, &kappa; = &tau;<sup>k</sup> / &sigma;<sub>applied</sub> , where &tau;<sup>k</sup> = tr(&rho;<sup>k</sup> &sigma;)
* Directional Strength to Stiffness Ratio (Y2E) = &epsilon;&#773; &tau;&#770; <sup>k</sup> / &tau; <sup>k</sup> , where &epsilon;&#773; is the nominal (average) strain and &tau;&#770; is the slip system strength of the most relevant slip system.



### Acknowledgements:

The Cornell University contributors' research was supported by the ONR under grant # N00014-16-1-3126, Dr. William Mullins Program Manager. The UC Santa Barbara contributors acknowledge the support of Office of Naval Research grants # N00014-19-1-2129 and N00014-18-1-2392 and National Science Foundation grant 1934641.
