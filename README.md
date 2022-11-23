# Density-Distribution-Plasma-Diagnostics

This project describes a numerical method used to construct the density distribution of a gas/plasma from the Mach Zehnder type interferometer fringes. Refractive index changes by gas molecules and plasma are characterized by the shifts of interfering fringes, which creates a 2D phase map. The gas/plasma density distribution is reconstructed from the phase map gradient, as described in the inverse Abel transformation method.

There are two major steps in implementing the scheme: 
1.	Extracting background fringes and the fringe shifts introduced by the gas/plasma, the difference between the gas/plasma fringe shifts and the background fringes are mapped into a two-dimensional phase distribution.
2.	Locating the symmetry axis from the peaks of the fringe; the gas/plasma density distribution is determined using the Inverse Abel Transformation from the phase gradient in the two-dimensional phase map.

The implementation details are shown in the following flowchart.
![image](https://user-images.githubusercontent.com/37091297/203460664-dbe0cb66-db83-48e2-bdbb-29aaac29309c.png)




Reference:
1.	Buckingham, A. D., & Graham, C. (1974). The density dependence of the refractivity of gases. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences, 337(1609), 275-291.

