# Node Explorer (3–10 nodes) extension

This package retains RVE Studio mechanics. It adds a discoverable custom geometry plugin called **NODESELECT** with a **NodeCount** parameter restricted to 3 through 10 inclusive. The eight seeded layouts correspond to an actual finite-candidate optimisation against the target `data/C_target_6x6_MPa.csv`. Note that some corner nodes are periodic images of the same junction.

In the GUI: Geometry Library -> activate **Node-count selector (3–10)**, set NodeCount minimum and maximum both to your chosen integer, and select StretchX, StretchY, Shear ranges. Use Generate DB and Large Batch Analysis as usual. The solver and Single RVE Designer can also load the `.rvejson` examples from the separate data package. This new selector is opt-in and does not change the solver or impose a global limit on manually imported projects.

The geometry plugin applies an affine transformation to the chosen fixed-connectivity base graph. For independent node-coordinate and wall-connectivity optimisation, use the reproducible research scripts `optimise_nodes.py`, `local_refinement.py`, and `refine_builtins.py` from the companion research package.

The target is a **minimum energy stiffness**, `CH-C_target` positive semidefinite, in Voigt [11,22,33,23,13,12]. This is NOT entrywise equality and many lightweight designs show a large matrix gap. Both cylindrical and free-moment results are checked. The stability check is a separate preliminary SSSS plate/Ritz screen for three illustrative applied-stress cases, not validation of structural buckling capacity.

To run the direct API example: `python examples/node_count_selector_demo.py` from the extracted project directory (after installing requirements).
