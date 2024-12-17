# Radar Simulations

This is a package for simulating the Spectral Polarimetric Variables of a Cloud Radar in moderate rain Conditions.
The spectral variables:

*  differential reflectivity (sZ<sub>DR</sub>)
*  differential phase (sδ<sub>hv</sub>)
*  co-polar correlation coefficient (sρ<sub>hv</sub>)

are generated taking into account the noise and the turbulence of real measurements. 
The **Documentation** and the **Theoretical Overview** are uploaded in pdf format and contain all the necessary information for the user. 

More information and content about the right fit of the Gamma drop size distribution and the second method of turbulence inclusion by using the covariance matrix, will be added after the publication of the paper *"Simulations of Spectral Polarimetric Variables measured in rain at W-band"* that was submitted to ACP on 14 Oct 2024.

## Licence

The project is available under the GNU Lesser General Public License v3.0. Please publish any changes you make, or better yet, contribute them to this repository so everyone can benefit. However, you can use this project as a library without publishing your source code.

The license is available in the `COPYING` and `COPYING.LESSER files`.

## References

-<span style="color:gray">Tsikoudi I., Battaglia A., Unal C., and Marinou E. *"Simulations of Spectral Polarimetric Variables measured in rain at W-band"*, Atmospheric Chemistry and Physics, submitted 2024. </span>

-This code is using the T-matrix scattering code from the module [pytmatrix](https://github.com/jleinonen/pytmatrix/wiki).

-Leinonen, Jussi. *"High-level interface to T-matrix scattering calculations: architecture, capabilities and limitations."* Optics express 22.2 (2014): 1655-1660.

-Mishchenko, Michael I., Joop W. Hovenier, and Larry D. Travis. *"Light scattering by nonspherical particles: theory, measurements, and applications."* Measurement Science and Technology 11.12 (2000): 1827-1827.


## Acknowledgement

This code was developed with the support by the Hellenic Foundation for Research and Innovation (H.F.R.I.) under the  “3rd Call for H.F.R.I. Research Projects to support Post-Doctoral Researchers” (Project Acronym: REVEAL, Project Number: 07222). 

## Contact

Ioanna Tsikoudi, National observatory of Athens (NOA) - jtsik@noa.gr
