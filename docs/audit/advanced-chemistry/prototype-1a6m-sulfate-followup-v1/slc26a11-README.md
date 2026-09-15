# slc26a11-ion-binding
## Calculating dissociation constants (Kd) from equilibrium MD simulations.

When calculating Kd from MD simulations, as with experimental methods, one must consider two entities, the ligand (chloride ion) and substrate (SLC26A11).

Usually, in experiments, Kd is calculated as the concentration of the ligand which which yields the half-maximal observable (e.g., florescence) when the substrate is in excess. The equilibrium constant Kd, which is defined as K_off/K_on, can easily be found, while the component rate-constants (K_off & K_on) are difficult to accquire.

In simulations, we often have a single substrate (1 protein) with an excess of ligand (e.g., chloride). The high spatial- and temporal resolution of MD simulations allows one to directly calculate the kinetic rate constants (K_off & K_on) directly from counting binding events, which is exactly what we do here.

Since we have a single molecule of substrate, we adopt a "protein-centric" view of ion-binding, wherein we don't care about the identity of the individual chloride, but only IF an chloride is bound or not. Such trajectories can be rather noisy, thus we utilize a Schmitt-Trigger to denoise the data, the boundaries of the Schmitt-Trigger are determined directly from the minimum distance trajectory. 

**ion_binding_analysis.py** was tested with MDAnalysis 2.5.0 (https://www.mdanalysis.org/), NumPy 1.26.4 (https://numpy.org/), and tqdm 4.67.1 (https://tqdm.github.io/).

The analyses are broken up into 3 steps:

	1. Generate a protein-centric minimum distance trajectory. At each frame the distance of the nearest
	   chloride to the anion binding pocket (ABP) is appended to the trajectory `ion_binding_analysis.pcBoundST()`.
	   Plotting these trajectories as histograms reveals a "minimum distance radial distribution".

	2. From the "minimum distance radial distribution", we identify a lower-bound (Lb) and an upper-bound (Ub)
	   for a Schmitt-Trigger, which is utilized in step-3. Lb would be the distance at which the first valley
	   following the "bound" peak appears. Ub would then be chosen as up-to 2 Å more than LB – specific
	   selections of Lb and Ub are system dependent.

	3. Finally, K_off, K_on, and Kd can be calculated using the minimum distance trajectory from step-1, and
	   Ub/Lb from step-2 `ion_binding_analysis.calkKd()`


**ion_binding_analysis.pcBoundST()**

Inputs:

	- uList		: Python List
				List containing MDAnalysis (MDA) Universe objects
	- ABPsel	: String
				MDA atomselection for the substrate-residues whose center-of-geometry is the binding-site
	- otherSel	: String
				MDA atomselection for the ligand. Default "name CL"
	- nSel		: Int
				Number of 'otherSel' considered for minimum-distance calculations. Default 1
	- lag 		: Int
				Number of frames to skip in-between distance evaluation

Outputs:

	- OUT		: Python List
				List containing minimum-distance trajectories. ion-binding-analysis.pcBoundST() assumes that
				the protein has two independent chains ('A', 'B'). Thus if the input list contains 1 universe,
				then `OUT` will contain two minimum-distance trajectories

**ion_binding_analysis.calcKd()**

Inputs:

	- trajList	: Python List
				The output list from ion_binding_analysis.pcBoundST(), i.e., `OUT`
	- dt 		: Float
				The time-step of the trajectory in nanoseconds. Default 0.1
	- lag		: Int
				Number of frames to skip in-between bound-state evaluation. Default 1
	- Lb		: Float
				Lower bound distance in Å, all values at or below Lb are considered 'Bound – 1'. Default 3
	- Ub		: Float
				Upper bound distance in Å, all values at or above Ub are considered 'Unbound – 0'. Default 5
	- BulkC		: Float
				Bulk concentration of the "otherSel" in Molar units. Default 0.2

Outputs:

	- Kds		: Numpy Array
				Array of [Kd, K_on, K_off] calculated for each minimum distance trajectory.
	- ONs		: Numpy Array
				Array of dwell-times for the Bound state.
	- OFFs		: Numpy Array
				Array of dwell-times for the Unbound state.
	- outTrajList	: Python List
				List of filtered state-trajectories, i.e., trajectory of 1's and 0's used for the Kd calculation.

Example Usage:

    import MDAnalysis as mda
	import ion_binding_analysis as iba
	import matplotlib.pyplot as plt

	# create a list of MDA universes – can also be a list of 1 universe.
	uList = []
	uList.append(mda.Universe("path/to/structure file", "path/to/coordinate file(s)"))

	# define your binding pocket and ligand
	ABPsel = "protein and resid XXX"
	otherSel = "name CL"

	# calculate minimum-distance trajectories
	trajList = iba.pcBoundST(uList=uList, ABPsel=ABPsel, otherSel=otherSel)

	# plot the minimum-distance radial distribution to identify Lb and Ub
	totalTraj = []
	for traj in trajList:
    	totalTraj.extend(traj[:,0])
	plt.hist(totalTraj)
	(example: Lb=3, Ub=5)

	# calculate Kd
	Kds, ONs, OFFs, outTrajList = iba.calcKd(trajList=trajList, dt=0.1, lag=1, Lb=3, Ub=5, BulkC=0.2)
	
# citation
If you make use of this data, please cite the corresponding paper:
```
@Article{Kuhn2026,
  author       = {Kuhn, Benedikt T. and Kovermann, Peter and Haddad, Bassam G. and Rasmussen, Tim and Hove, Tamsanqa T. and Bungert-Plümke, Stefanie and Böttcher, Bettina and Machtens, Jan-Philipp and Fahlke, Christoph and Geertsma, Eric R.},
  date         = {2026-07},
  journaltitle = {Nature Communications},
  title        = {SLC26A11 is an atypical solute carrier with dual transport-channel function mediating lysosomal sulfate transport},
  doi          = {10.1038/s41467-026-75749-4},
  issn         = {2041-1723},
  number       = {1},
  volume       = {17},
  publisher    = {Springer Science and Business Media LLC},
}
```
