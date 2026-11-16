# Archived FE documents

Superseded by [../FE_VERIFICATION.md](../FE_VERIFICATION.md) (RVE Studio 1.11). They describe front ends that
were removed from the GUI and the command line in 1.11 (the single-run Abaqus page, `fe-dbc`, `fe-campaign`,
`fe-study`, `verify`, other periodic domains, 21-case energy recovery). The Python functions they mention
remain in `rvestudio.abaqus_thesis_dbc`, `rvestudio.abaqus_homogenization` and `rvestudio.fe_validation`.
Their top/bottom conditions are those of 1.10: since 1.11 the DBC model adds four cases with U3 = 0 on the
top and bottom faces (3D terms) and the PBC model has no E33 state (U3+ - U3- = 0).
