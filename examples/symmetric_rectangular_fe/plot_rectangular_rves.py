"""Plot the exact geometric rectangular FE windows (not FE results)."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from rvestudio.io import load_project
from rvestudio.rectangular_selection import inspect_rectangular_rve
from rvestudio.abaqus_homogenization import build_rectangular_supercell

def make_plot(out):
    base=Path(__file__).resolve().parent
    paths=[base/'SQUARE_symmetric_0_90_0.rvejson',base/'MIXED_symmetric_0_90_0.rvejson',base.parent/'fe_validation'/'P00050.rvejson']
    titles=['SQUARE: symmetric rectangular','MIXED: symmetric rectangular','P00050: unrestricted rectangular']
    fig,axes=plt.subplots(1,3,figsize=(13.8,5.1),layout='constrained')
    for ax,path,title in zip(axes,paths,titles):
        primitive,mats,stacks=load_project(path)
        if 'P00050' in path.name:
            m,d=build_rectangular_supercell(primitive,repeats_x=1,repeats_y=1)
            diag='periodic, not mirror-certified'
        else:
            m,d=inspect_rectangular_rve(primitive,stacks)
            diag='mirror-checked in x and y'
        B=np.column_stack((m.a1,m.a2))
        uv=np.linalg.solve(B,np.asarray(m.nodes).T).T
        for idx,e in enumerate(m.elements):
            p,q=uv[e.node_i],uv[e.node_j]
            ax.plot([p[0],q[0]],[p[1],q[1]],linewidth=2.3)
            ax.scatter([p[0],q[0]],[p[1],q[1]],s=12)
        bd=np.array([[0,0],[1,0],[1,1],[0,1],[0,0]])
        ax.plot(bd[:,0],bd[:,1],color='black',linewidth=1.7)
        ax.set_xlim(-.1,1.1);ax.set_ylim(-.12,1.12);ax.set_aspect('equal');ax.grid(alpha=.16)
        ax.set_title(title,fontsize=11,weight='bold')
        ax.set_xlabel(f"x / {np.linalg.norm(m.a1):.2f} mm")
        ax.set_ylabel(f"y / {np.linalg.norm(m.a2):.2f} mm")
        ax.text(.5,-.18,diag,transform=ax.transAxes,ha='center',fontsize=9)
    fig.suptitle('Commensurate rectangular periodic FE domains (wall midsurfaces, top view)\nAll S4 shells extend vertically through the modelled core height; no Abaqus FE solution shown',fontsize=12)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=165,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':
    make_plot(Path(__file__).resolve().parent/'rectangular_fe_geometries.png')
