# RVE Studio — Analysis and buckling scope

The linear wall-buckling workflow is retained from the preceding application release. Each wall's in-plane membrane resultants `[Nbar11,Nbar33,Nbar13]` are recovered from the homogenisation solution and used in a simply-supported laminated-plate eigenvalue screening calculation. The whole combined load scales with the eigenfactor. Unsymmetric laminates use the documented Schur-reduced bending matrix `D* = D − B A^-1 B` by default.

The referenced article by Sriharan et al. (2024), Engineering Structures 299, 117106, provides the **homogenisation reference**. It is not presented as the source of the separate wall-buckling screening approximation.

Publication Eq. (5) defines the reduced shell transverse-shear order `[13,23]`, giving Q44=C1313 and Q55=C2323. Internally the 3D tensor retains Voigt `[11,22,33,23,13,12]`. This distinction is essential when interpreting databases from earlier releases.

With the default Euler–Bernoulli walls the periodic-RVE implementation reproduces all printed paper stiffness values within 0.3%; the Timoshenko wall option lowers Q66 by up to 14.5% relative to the paper. Refer to `VALIDATION_REPORT.txt` and `data/benchmark_all14_detail.csv` for complete values and avoid interpreting a configurable benchmark threshold as independent FE validation.
