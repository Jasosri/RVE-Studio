"""Plot *unanalysed* Abaqus S4 shell mesh geometry, not FE stress results."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from rvestudio.io import load_project
from rvestudio.abaqus_homogenization import resolve_periodic_domain, PeriodicHomogenizationSettings
from rvestudio.abaqus import _build_shell_mesh


def render(project, output, mesh_size=10.0):
    model, _, _ = load_project(project)
    mod, _ = resolve_periodic_domain(model, PeriodicHomogenizationSettings(
        mesh_size=mesh_size,domain_mode='lattice_supercell',
        supercell_repeat_x=2,supercell_repeat_y=2))
    nodes, wall_elements, _, _ = _build_shell_mesh(mod, mesh_size)
    faces=[]
    for elems in wall_elements.values():
        for _, conn in elems:
            faces.append([nodes[n] for n in conn])
    fig=plt.figure(figsize=(8.2,6.5));ax=fig.add_subplot(111,projection='3d')
    group=Poly3DCollection(faces,alpha=.7,edgecolors=(.14,.25,.38,.22),linewidths=.18)
    group.set_facecolor((.38,.66,.86,.42));ax.add_collection3d(group)
    xyz=list(nodes.values())
    ax.auto_scale_xyz([p[0] for p in xyz],[p[1] for p in xyz],[p[2] for p in xyz])
    ax.set_box_aspect((1,1,1.1));ax.view_init(elev=21,azim=-58)
    ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)');ax.set_zlabel('z (mm)')
    ax.set_title(f'{Path(project).stem}: 2×2 periodic S4 shell mesh, h={mesh_size:g} mm\nFE deck preparation; no stress or displacement solution')
    fig.tight_layout();fig.savefig(output,dpi=155);plt.close(fig)
    return len(faces)

if __name__=='__main__':
    root=Path(__file__).resolve().parents[0]
    import sys
    for label in ('P00050','P00148','P00004'):
        out=Path(sys.argv[1]) if len(sys.argv)>1 else root
        out.mkdir(exist_ok=True,parents=True)
        print(label, render(root / f'{label}.rvejson',out / f'{label}_S4_mesh_2x2.png'))
