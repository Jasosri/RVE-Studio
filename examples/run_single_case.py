from rvestudio.examples_data import paper_hexagonal_case
from rvestudio.homogenization import homogenize_rve

model,materials,laminates=paper_hexagonal_case(1)
out=homogenize_rve(model,materials,laminates)
print('C^H [MPa]')
print(out.CH)
print('\nQplane')
print(out.props.Qplane)
print('\nQshear')
print(out.props.Qshear)
print('\nRelative density:',out.relative_density)
print('Mass [kg]:',out.mass_kg)
