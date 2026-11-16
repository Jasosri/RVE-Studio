from pathlib import Path
import sys,zipfile,io,math,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,LogNorm
from matplotlib.ticker import MaxNLocator
P=Path(__file__).resolve().parent; F=P/'figures'; F.mkdir(exist_ok=True)
for candidate in (P/'RVE_Studio_Node_Explorer',P.parent.parent,P.parent.parent.parent):
 if (candidate/'rvestudio').is_dir(): sys.path.insert(0,str(candidate)); break
from rvestudio.geometry import BUILTIN_GEOMETRIES
import run_cylindrical_study as study
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'figure.dpi':150,'savefig.dpi':170,'axes.spines.top':False,'axes.spines.right':False})
BLUE='#1F497D';CYAN='#2BA6A0'; ORANGE='#D18A32';RED='#B94D45';GREY='#687684'; PALETTE={'HEX':BLUE,'SQUARE':GREY,'MIXED':ORANGE,'TRI':CYAN}

def save(fig,name):
 fig.savefig(F/name,bbox_inches='tight',facecolor='white');plt.close(fig)
def load(n):return pd.read_csv(P/'prior_cyl'/n)
atlas=load('atlas_solved_paired_2400.csv');fits=load('controlled_scaling_fits.csv');robust=load('robust_optimisation_all_designs.csv');loads=load('load_interaction_critical_walls.csv');custom=load('custom_matched_effect.csv');custompairs=load('custom_vs_square_matched_120.csv')
bench=load('14config_paper_cyl.csv')
initial=pd.read_csv(P/'results/initial_screen.csv')
refined=pd.read_csv(P/'results/stiffness_refined.csv')
bk=pd.read_csv(P/'results/buckling_refined.csv');feas=bk[bk.Status=='PASS'].copy()
T=np.loadtxt(P/'results/C_target_MPa.csv',delimiter=',');C=np.loadtxt(P/'results/C_achieved_P00050_MPa.csv',delimiter=',');Ch=np.loadtxt(P/'results/C_achieved_P00148_MPa.csv',delimiter=',')

# F01: schematics
fig,axs=plt.subplots(2,2,figsize=(10,6.6));
for ax,gid in zip(axs.flat,study.COLORS):
 pr={'HEX':(52.5,52.5),'SQUARE':(52.5,52.5),'MIXED':(45,2.),'TRI':(52.5,1.)}[gid];geo=study.geometry_meta(gid,pr)
 for i,(a,b) in enumerate(geo.connectivity):
  xs=geo.nodes[[a,b],0];ys=geo.nodes[[a,b],1];ax.plot(xs,ys,color=PALETTE[gid],linewidth=3,alpha=.9,solid_capstyle='round');ax.text(xs.mean(),ys.mean(),str(i+1),fontsize=8,color='#394D5B')
 ax.scatter(geo.nodes[:,0],geo.nodes[:,1],c=BLUE,s=23,zorder=2)
 for k,xy in enumerate(geo.nodes):ax.annotate(str(k+1),xy,xytext=(4,3),textcoords='offset points',fontsize=7)
 ax.set_aspect('equal');ax.set_title(f'{gid}  |  {len(geo.nodes)} nodes · {len(geo.connectivity)} walls',weight='bold');ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]');ax.grid(alpha=.14)
fig.suptitle('Four built-in graph families (geometric examples)',fontsize=13,fontweight='bold');fig.tight_layout();save(fig,'F01_builtin_geometries.png')
# F02 benchmark
b=bench.copy();b['pct']=100*(b.PaperCyl_MPa-b.Published_MPa)/b.Published_MPa
fig,axs=plt.subplots(1,2,figsize=(10.4,3.6))
for prop,grp in b.groupby('Property',sort=False): axs[0].plot(grp.Config,grp.pct,'o-',markersize=3,label=prop)
axs[0].axhline(1,ls='--',color=GREY,lw=.8);axs[0].axhline(-1,ls='--',color=GREY,lw=.8)
axs[0].set_xlabel('Published configuration');axs[0].set_ylabel('Signed difference [%]');axs[0].set_title('All 14 configurations · six shell terms');axs[0].legend(ncol=3,fontsize=7);axs[0].grid(alpha=.18)
ss=b.groupby('Property').pct.apply(lambda x:np.max(np.abs(x))).reindex(['Q11','Q22','Q12','Q66','Q44','Q55']);axs[1].bar(ss.index,ss.values,color=[RED if x=='Q66' else BLUE for x in ss.index]);axs[1].set_ylabel('Maximum absolute difference [%]');axs[1].set_title('Unresolved Q66 benchmark discrepancy');axs[1].grid(axis='y',alpha=.15)
fig.tight_layout();save(fig,'F02_journal_benchmark_cyl.png')
# F03 atlas
fig,axs=plt.subplots(1,2,figsize=(10.4,4))
for gid,g in atlas.groupby('GeometryID'):
 axs[0].scatter(g.AverageDensityKgM3,g.Qmin_paper,s=7,alpha=.22,color=PALETTE[gid],label=gid)
axs[0].set_xlabel('Effective density [kg/m³]');axs[0].set_ylabel('min(Q11, Q22) [MPa]');axs[0].set_title('2,400-design cylinder-only atlas');axs[0].legend();axs[0].grid(alpha=.12)
median=atlas.groupby('GeometryID').SpecificQmin_paper.median().reindex(study.COLORS)
axs[1].bar(median.index,median.values,color=[PALETTE[g] for g in median.index]);axs[1].set_ylabel('Median Qmin / effective density');axs[1].set_title('Specific-property comparison · unpaired');axs[1].grid(axis='y',alpha=.12)
fig.tight_layout();save(fig,'F03_atlas_specific_stiffness.png')
# F04 scaling
ff=fits[(fits.Mode=='paper')&(fits.Metric=='Q66_MPa')].copy();fqq=fits[(fits.Mode=='paper')&(fits.Metric=='Qmin_MPa')].copy()
fig,axs=plt.subplots(1,2,figsize=(10.2,3.7))
for ax,dat,name in [(axs[0],ff,'Q66'),(axs[1],fqq,'Qmin')]:
 means=dat.groupby('GeometryID').n.agg(['median','min','max']).reindex(study.COLORS)
 yerr=np.vstack([means['median']-means['min'],means['max']-means['median']]);ax.bar(means.index,means['median'],yerr=yerr,capsize=4,color=[PALETTE[g] for g in means.index]);ax.set_ylim(bottom=0);ax.set_title(f'Power-law exponent in {name} ~ (ρ*)ⁿ');ax.set_ylabel('Exponent n');ax.grid(axis='y',alpha=.15)
fig.tight_layout();save(fig,'F04_density_scaling.png')
# F05 cylindrical robust screening - prior
col_nom='NominalFeasible_paper';col_cons='ConservativeFeasible_paper'
print('ROBUST COLUMNS', [x for x in robust.columns if 'Feasible' in x or 'Pareto' in x][:20])
fig,axs=plt.subplots(1,2,figsize=(10.4,3.9))
for gid,g in robust.groupby('GeometryID'):
 axs[0].scatter(g.AverageDensityKgM3,g.Qmin_paper,s=6,alpha=.18,color=PALETTE[gid],label=gid)
axs[0].set_xlabel('Density [kg/m³]');axs[0].set_ylabel('Qmin [MPa]');axs[0].set_title('Nominal performance and design trade-off');axs[0].grid(alpha=.15)
rr=robust[['GeometryID','Q11_paper','Q22_paper']].copy();rr['ratio']=rr.Q11_paper/rr.Q22_paper
for gid,g in rr.groupby('GeometryID'):
 axs[1].hist(g.ratio,bins=np.linspace(.35,4,45),alpha=.45,label=gid,color=PALETTE[gid])
axs[1].axvline(1.5,color=RED,ls='--',label='Inverse target 1.5');axs[1].set_xlim(.35,4);axs[1].set_xlabel('Q11/Q22');axs[1].set_ylabel('Design count');axs[1].set_title('Inverse target-proportion screening');axs[1].legend(fontsize=7)
fig.tight_layout();save(fig,'F05_robust_inverse_nominal.png')
# F06 critical wall, cylindrical loads
lc=loads[loads.Mode=='paper'].copy();fig,axs=plt.subplots(1,2,figsize=(10.4,4.2));
for gid,g in lc.groupby('GeometryID'):
 axs[0].scatter(g.Sigma11_MPa,g.BucklingMultiplier,s=22,alpha=.6,label=gid,color=PALETTE.get(gid,GREY))
axs[0].axhline(1.5,color=RED,ls='--',lw=1,label='Minimum factor');axs[0].set_yscale('log');axs[0].set_xlabel('Applied σ11 [MPa]');axs[0].set_ylabel('Buckling load factor (log)');axs[0].set_title('Load-dependent wall-stability screen');axs[0].legend(ncol=2,fontsize=7);axs[0].grid(alpha=.13)
ct=lc.groupby(['GeometryID','CriticalWall']).size().unstack(fill_value=0).reindex(study.COLORS).fillna(0);ct.plot(kind='bar',stacked=True,ax=axs[1],colormap='tab20c',legend=False);axs[1].set_ylabel('Number of governing-wall observations');axs[1].set_xlabel('Architecture');axs[1].set_title('Critical wall switches with loading');axs[1].grid(axis='y',alpha=.14)
fig.tight_layout();save(fig,'F06_buckling_load_interaction.png')
# F07 custom matched, cylindrical only
cm=custom[custom.Mode=='paper'].copy();fig,axs=plt.subplots(1,2,figsize=(10.2,4.))
for ax,xx,yy,ttl in [(axs[0],'Qmin_MPa_SQUARE','Qmin_MPa_BRACED','Weak-direction stiffness'),(axs[1],'Q66_MPa_SQUARE','Q66_MPa_BRACED','In-plane shear stiffness')]:
 ax.scatter(cm[xx],cm[yy],s=13,alpha=.5,color=CYAN);lims=[min(cm[xx].min(),cm[yy].min()),max(cm[xx].max(),cm[yy].max())];ax.plot(lims,lims,'--',color=RED,lw=1);ax.set_xlim(lims);ax.set_ylim(lims);ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Square control [MPa]');ax.set_ylabel('Diagonal braced [MPa]');ax.set_title(ttl);ax.grid(alpha=.13)
fig.tight_layout();save(fig,'F07_custom_matched_cyl.png')
# F08 optimisation funnel
fig,axs=plt.subplots(1,2,figsize=(10.4,4.))
counts=[len(initial),int((initial.Adequacy>=1.05).sum()),len(refined),len(bk),len(feas)]
labels=['Initial\n4,128','At ρ*=0.15:\nstiffness','Thickness\nrefinements','Three-load\nshortlist','Passing all\nconstraints'];bars=axs[0].bar(range(len(counts)),counts,color=[BLUE,CYAN,ORANGE,GREY,BLUE]);axs[0].set_xticks(range(len(counts)),labels);axs[0].set_ylabel('Number of candidates');axs[0].set_title('Finite-sample search funnel');
for bar,v in zip(bars,counts):axs[0].text(bar.get_x()+bar.get_width()/2,v+25,str(v),ha='center',fontsize=8)
for gid,g in initial[initial.Status=='OK'].groupby('Geometry'):
 axs[1].scatter(g.Density,g.Adequacy,s=4,alpha=.2,label=gid,color=PALETTE[gid]);
axs[1].axhline(1.05,ls='--',color=RED);axs[1].set_xlabel('Density at ρ*=0.15 [kg/m³]');axs[1].set_ylabel('Minimum generalised stiffness factor');axs[1].set_yscale('log');axs[1].set_title('Initial target-C screening');axs[1].legend(fontsize=7)
fig.tight_layout();save(fig,'F08_search_funnel.png')
# F09 best by geometry and ply
fig,axs=plt.subplots(1,2,figsize=(10.4,3.85));gb=feas.groupby('Geometry').Density.min().reindex(study.COLORS);axs[0].bar(gb.index,gb.values,color=[PALETTE[g] for g in gb.index]);axs[0].set_ylabel('Minimum feasible density [kg/m³]');axs[0].set_title('Best in screened subset by geometry');axs[0].grid(axis='y',alpha=.15)
piv=feas.pivot_table(index='MaterialScheme',columns='PlyCount',values='Density',aggfunc='min').reindex(['Glass','Boron','Hybrid_GBG','Hybrid_BGB'])
im=axs[1].imshow(piv.values,cmap='YlGnBu_r',aspect='auto',vmin=np.nanmin(piv.values),vmax=np.nanmax(piv.values));axs[1].set_xticks(np.arange(len(piv.columns)),piv.columns);axs[1].set_yticks(np.arange(len(piv.index)),piv.index);axs[1].set_xlabel('Ply count');axs[1].set_title('Lightest feasible density [kg/m³]');
for i in range(len(piv.index)):
 for j in range(len(piv.columns)):
  x=piv.iat[i,j];axs[1].text(j,i,('-' if np.isnan(x) else f'{x:.1f}'),ha='center',va='center',fontsize=9,color='black')
fig.colorbar(im,ax=axs[1],shrink=.9);fig.tight_layout();save(fig,'F09_geometry_material_ply_choice.png')
# F10 feasibility tradeoff
fig,axs=plt.subplots(1,2,figsize=(10.4,4.));
for gid,g in bk[bk.Status=='PASS'].groupby('Geometry'):
 axs[0].scatter(g.Density,g.MatrixGap,s=22,alpha=.75,label=gid,color=PALETTE[gid])
axs[0].scatter([feas.iloc[0].Density],[feas.iloc[0].MatrixGap],color=RED,marker='*',s=160,zorder=4,label='Minimum mass');axs[0].set_xlabel('Density [kg/m³]');axs[0].set_ylabel('Relative Frobenius C-matrix gap');axs[0].set_title('Matching target vs satisfying lower bound');axs[0].legend(fontsize=7);axs[0].grid(alpha=.15)
for mat,g in feas.groupby('MaterialScheme'):
 axs[1].scatter(g.ArealMass,g.TotalCostPerM2,s=23,alpha=.75,label=mat)
axs[1].set_xlabel('Areal mass [kg/m²]');axs[1].set_ylabel('Illustrative material cost index / m²');axs[1].set_title('Mass–cost screening');axs[1].legend(fontsize=7);axs[1].grid(alpha=.15)
fig.tight_layout();save(fig,'F10_mass_cost_matrix_matching.png')
# F11 target/achieved matrices
fig,axs=plt.subplots(1,3,figsize=(12,4.0)); vmax=max(T.max(),C.max());
for ax,mat,ttl in zip(axs,[T,C,Ch],['Required target C','Minimum-mass P00050','Lower-cost hybrid P00148']):
 im=ax.imshow(mat,cmap='Blues',vmin=0,vmax=vmax)
 for i in range(6):
  for j in range(6):
   if abs(mat[i,j])>=1.0:ax.text(j,i,f'{mat[i,j]:.0f}',ha='center',va='center',fontsize=7,color='white' if mat[i,j]>vmax*.5 else '#14314a')
 ax.set_xticks(range(6),['11','22','33','23','13','12']);ax.set_yticks(range(6),['11','22','33','23','13','12']);ax.set_title(ttl)
fig.colorbar(im,ax=axs,fraction=.03,pad=.02,label='Stiffness [MPa]');save(fig,'F11_target_and_achieved_matrices.png')
# F12 selected geometry + laminate stack
fig,axs=plt.subplots(1,2,figsize=(10.4,4.));m=next(z for z in study.allmeta if z['ID']=='P00050');g=study.geometry_meta('HEX',(m['Param1'],m['Param2']));ax=axs[0]
for j,(a,b) in enumerate(g.connectivity):
 ax.plot(g.nodes[[a,b],0],g.nodes[[a,b],1],lw=6,color=ORANGE,solid_capstyle='round');xy=g.nodes[[a,b]].mean(axis=0);ax.text(xy[0],xy[1],f'wall {j+1}',fontsize=9,ha='center',va='center',bbox={'facecolor':'white','edgecolor':'none','alpha':.86})
ax.scatter(g.nodes[:,0],g.nodes[:,1],c=BLUE,s=28);ax.set_aspect('equal');ax.grid(alpha=.18);ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]');ax.set_title('Optimal sampled HEX unit graph');
ax=axs[1]
for i,a in enumerate([0,45,0]):
 ax.barh(i,1,left=0,height=.86,color=[BLUE,CYAN,ORANGE][i]);ax.text(.5,i,f'Ply {i+1}: {a}°  Boron/Epoxy',ha='center',va='center',color='white',fontsize=10,weight='bold')
ax.invert_yaxis();ax.set_yticks([]);ax.set_xticks([]);ax.set_xlim(0,1);ax.set_title('Optimised 3-ply symmetric [0/45/0]');ax.text(.5,3.0,'Equal-ply thickness; wall thickness 1.816 mm',ha='center',fontsize=9)
fig.tight_layout();save(fig,'F12_selected_geometry_layup.png')
# F13 three loads, critical wall
fig,axs=plt.subplots(1,2,figsize=(10.4,3.6)); ids=['P00050','P00148','P03416','P02653']; cols=[BLUE,CYAN,ORANGE,GREY];xx=np.arange(3);wd=.19
for i,cid in enumerate(ids):
 row=bk[bk.ID==cid].iloc[0];v=[row[f'Buck_{z}'] for z in study.LOADS];axs[0].bar(xx+(i-1.5)*wd,v,wd,label=cid,color=cols[i])
axs[0].axhline(1.5,ls='--',color=RED,label='Criterion');axs[0].set_xticks(xx,['Uniaxial','Biaxial','Combined']);axs[0].set_ylabel('Linear buckling multiplier');axs[0].set_title('Selected designs under three load states');axs[0].legend(fontsize=7,ncol=3);axs[0].grid(axis='y',alpha=.13)
for i,cid in enumerate(ids):
 row=bk[bk.ID==cid].iloc[0];v=[row[f'Wall_{z}'] for z in study.LOADS];axs[1].plot(xx,v,'o-',label=cid,color=cols[i])
axs[1].set_xticks(xx,['Uniaxial','Biaxial','Combined']);axs[1].set_ylabel('Governing wall ID');axs[1].set_title('Critical-wall switching');axs[1].yaxis.set_major_locator(MaxNLocator(integer=True));axs[1].grid(alpha=.16);axs[1].legend(fontsize=7)
fig.tight_layout();save(fig,'F13_buckling_selected.png')
# F14 constraint margin bar
f=feas.sort_values('Density').iloc[0];values=[f.Adequacy/1.05,f.Buck_uniaxial/1.5,f.Buck_biaxial/1.5,f.Buck_combined/1.5,50/f.TotalCostPerM2,f.MinThickness/.35,f.MinPlyThickness/.08,1.9/f.ThicknessRatio]
labels=['Stiffness','Buckling U','Buckling B','Buckling C','Cost cap','Wall t min','Ply t min','t max/min'];fig,ax=plt.subplots(figsize=(9.5,3.8));bars=ax.barh(labels,values,color=[ORANGE if v<1.1 else BLUE for v in values]);ax.axvline(1,color=RED,ls='--',lw=1.1);ax.invert_yaxis();ax.set_xlabel('Achieved / required feasibility ratio (>=1 passes)');ax.set_title('Constraint margins for minimum-mass P00050');ax.grid(axis='x',alpha=.16)
for i,v in enumerate(values):ax.text(v+.08,i,f'{v:.2f}',va='center',fontsize=8)
fig.tight_layout();save(fig,'F14_constraint_margins.png')
# F15 custom workflow
fig,ax=plt.subplots(figsize=(10.4,1.9));ax.axis('off');steps=['Choose built-in\ngeometry','Choose material\n& plies','Construct periodic\nwall graph','Cylindrical\nhomogenisation','C-target +\nbuckling checks','Rank minimum\nmass solutions'];loc=np.linspace(.1,.9,len(steps))
for i,(x,txt) in enumerate(zip(loc,steps)):
 ax.text(x,.55,txt,ha='center',va='center',fontsize=8.5,color='white',weight='bold',bbox=dict(boxstyle='round,pad=.5',fc=BLUE if i%2==0 else CYAN,ec='none'))
 if i<len(steps)-1:ax.annotate('',xy=(loc[i+1]-.068,.55),xytext=(x+.068,.55),arrowprops=dict(arrowstyle='->',color=GREY,lw=1.7))
ax.set_title('Reproducible cylindrical-only optimisation workflow',weight='bold');save(fig,'F15_workflow.png')
print('FIGURES',len(list(F.glob('*.png'))),'benchmark within1',sum(abs(bench.PaperCyl_MPa/bench.Published_MPa-1)*100<=1),'n benchmark',len(bench),'median atlas',atlas.groupby('GeometryID').Qmin_paper.median().to_dict())
