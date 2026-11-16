"""RVE Studio — laminated composite cellular-core homogenization."""
__version__='1.0.0'
from .models import Material,Ply,Laminate,Element,RVEModel
from .homogenization import homogenize_rve
from .laminate import analyze_laminate,laminate_section,material_invariants,lamination_parameters
from .examples_data import paper_hexagonal_case,thesis_hexagonal_case,default_materials

from .load_recovery import recover_wall_loads, WallLoadRecoveryResult

from .buckling import laminated_plate_buckling, evaluate_wall_buckling

from .api import solve, solve_example, solve_project, summary, benchmark

# Abaqus input preparation and strict rectangular FE geometry selection.
from .abaqus_homogenization import PeriodicHomogenizationSettings, PBC_DOMAIN_SYMMETRIC_RECTANGULAR
from .rectangular_selection import inspect_rectangular_rve
