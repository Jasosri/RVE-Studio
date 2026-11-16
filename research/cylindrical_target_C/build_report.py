from pathlib import Path
import sys,io,zipfile,json,math
import numpy as np,pandas as pd
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH,WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_LINE_SPACING
from docx.enum.style import WD_STYLE_TYPE
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.section import WD_SECTION_START
P=Path(__file__).resolve().parent;OUT=P/'RVE_Studio_Cylindrical_Only_Five_Case_Studies_Manuscript.docx'
R=P/'results';F=P/'figures';initial=pd.read_csv(R/'initial_screen.csv');refined=pd.read_csv(R/'stiffness_refined.csv');bk=pd.read_csv(R/'buckling_refined.csv');ok=bk[bk.Status=='PASS'].copy();best=ok.nsmallest(1,'Density').iloc[0]
def load(n):return pd.read_csv(P/'prior_cyl'/n)
atlas=load('atlas_solved_paired_2400.csv');fits=load('controlled_scaling_fits.csv');rb=load('robust_optimisation_all_designs.csv');critical=load('load_interaction_critical_walls.csv');custom=load('custom_matched_effect.csv');rinv=pd.DataFrame()
bench=load('14config_paper_cyl.csv')
T=np.loadtxt(R/'C_target_MPa.csv',delimiter=',');C=np.loadtxt(R/'C_achieved_P00050_MPa.csv',delimiter=',');H=np.loadtxt(R/'C_achieved_P00148_MPa.csv',delimiter=',')

DOC=Document();sec=DOC.sections[0];sec.top_margin=Inches(.74);sec.bottom_margin=Inches(.7);sec.left_margin=Inches(.78);sec.right_margin=Inches(.78);sec.header_distance=Inches(.28);sec.footer_distance=Inches(.32)
styles=DOC.styles
styles['Normal'].font.name='Aptos';styles['Normal'].font.size=Pt(9.5);styles['Normal'].paragraph_format.space_after=Pt(6);styles['Normal'].paragraph_format.line_spacing=1.18
for n,sz,col in [('Title',20,(17,47,76)),('Heading 1',15,(20,57,91)),('Heading 2',11.8,(20,57,91)),('Heading 3',10,(18,72,104))]:
 st=styles[n];st.font.name='Aptos Display' if n!='Heading 3' else 'Aptos';st.font.size=Pt(sz);st.font.color.rgb=RGBColor(*col);st.font.bold=True;st.paragraph_format.space_before=Pt(13 if n!='Title' else 0);st.paragraph_format.space_after=Pt(6)
for style_name in ['Caption']:
 s=styles[style_name];s.font.name='Aptos';s.font.size=Pt(8);s.font.italic=True;s.font.color.rgb=RGBColor(74,84,95);s.paragraph_format.space_after=Pt(9)
# custom section lead
if 'Lead' not in styles:
 s=styles.add_style('Lead',WD_STYLE_TYPE.PARAGRAPH);s.base_style=styles['Normal'];s.font.size=Pt(10.3);s.font.color.rgb=RGBColor(50,68,80)

head=sec.header.paragraphs[0];head.text='RVE STUDIO v1   |   CYLINDRICAL-BENDING STUDIES';head.style='Caption';head.alignment=WD_ALIGN_PARAGRAPH.RIGHT
foot=sec.footer.paragraphs[0];foot.alignment=WD_ALIGN_PARAGRAPH.CENTER
foot.add_run('Computational research manuscript · 30 October 2026     |     ')
fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');foot._p.append(fld)

def para(t='',style=None,keep=False):
 p=DOC.add_paragraph(style=style or 'Normal');p.add_run(t)
 if keep:p.paragraph_format.keep_with_next=True
 return p

def title(t):DOC.add_heading(t,1)
def sub(t):DOC.add_heading(t,2)
def eq(t):
 p=DOC.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_before=Pt(4);p.paragraph_format.space_after=Pt(10);r=p.add_run(t);r.bold=True;r.font.name='Cambria Math';r.font.size=Pt(10.3)
 return p

def bullet(t):DOC.add_paragraph(t,style='List Bullet')

def caption(t):p=DOC.add_paragraph(style='Caption');p.add_run(t);return p

def fig(fname,num,txt,width=6.68):
 p=DOC.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER;p.paragraph_format.space_before=Pt(9);p.paragraph_format.keep_with_next=True
 p.add_run().add_picture(str(F/fname),width=Inches(width));caption(f'Figure {num}. {txt}')

def shade(c,fill):
 tcPr=c._tc.get_or_add_tcPr();x=OxmlElement('w:shd');x.set(qn('w:fill'),fill);tcPr.append(x)

def table(cols,rows,caption_text=None,widths=None):
 if caption_text:
  p=DOC.add_paragraph(style='Caption');p.add_run(caption_text);p.paragraph_format.keep_with_next=True
 t=DOC.add_table(rows=1,cols=len(cols));t.style='Table Grid';t.alignment=WD_TABLE_ALIGNMENT.CENTER
 for j,c in enumerate(cols):
  cell=t.rows[0].cells[j];cell.text=str(c);shade(cell,'1F497D')
  for pa in cell.paragraphs:
   for r in pa.runs:r.font.color.rgb=RGBColor(255,255,255);r.bold=True;r.font.size=Pt(8)
 for i,row in enumerate(rows):
  cells=t.add_row().cells
  for j,v in enumerate(row):
   cells[j].text=str(v);cells[j].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
   for pa in cells[j].paragraphs:
    pa.paragraph_format.space_after=Pt(2);pa.paragraph_format.space_before=Pt(2)
    for r in pa.runs:r.font.size=Pt(7.6)
   if i%2:shade(cells[j],'F1F5F9')
 for row in t.rows:
  for j,cell in enumerate(row.cells):
   if widths is not None:cell.width=Inches(widths[j])
 return t

def page():pass

p=DOC.add_paragraph(style='Title');p.text='High-Throughput Homogenisation and Minimum-Mass Inverse Design of Laminated Composite Cellular Cores';p.paragraph_format.space_after=Pt(9)
p=DOC.add_paragraph();r=p.add_run('RVE Studio v1 — cylindrical-bending formulation and five reproducible computational case studies');r.font.size=Pt(13);r.font.bold=True;r.font.color.rgb=RGBColor(36,91,128)
para('Jasotharan Sriharan',style='Lead')
para('Institute for Infrastructure and Environment, School of Engineering, The University of Edinburgh, Edinburgh, United Kingdom')
para('Research-software and computational mechanics methodology · Working manuscript · 30 October 2026')
DOC.add_paragraph('ABSTRACT',style='Heading 2')
para('This manuscript presents a cylindrical-bending-only computational investigation using RVE Studio v1, a Python implementation of strain-energy homogenisation for laminated composite cellular cores. Unlike previous experimental research drafts that contrasted cylindrical and free-moment strip assumptions, every homogenised stiffness and load-recovery result considered here uses the paper-style cylindrical laminated Timoshenko wall model. Five worked studies illustrate a 2,400-design topology–material atlas, robust multi-objective inverse design, combined-loading wall-buckling screens, custom periodic graph comparison, and a new target-C minimum-mass search over 4,128 combinations of built-in geometry parameters, material arrangements, ply counts (3, 5, 7, 9), and symmetric stacking sequences. The final search adopts a complete positive-definite 6×6 target stiffness matrix, a minimum 5% generalised strain-energy stiffness margin, explicit wall/ply manufacturing bounds, an illustrative material-cost ceiling, and linear wall-buckling factors of at least 1.5 under three prescribed loading states. Among 136 stability-screened candidates, 110 satisfy all constraints; the lightest sampled solution is a three-ply [0/45/0] boron/epoxy hexagonal core at 83.09 kg/m³ (12.46 kg/m² at 150 mm height), with minimum generalised stiffness factor 1.0505 and minimum screened buckling factor 3.686. The alternative five-ply hybrid design is only 0.55% heavier and has lower illustrative cost. The target tensor is exceeded in the strain-energy sense rather than replicated entry-by-entry, and the Q66 publication discrepancy and simplified buckling assumptions are explicitly retained as limitations.')
para('Keywords: composite cellular cores; periodic homogenisation; cylindrical bending; target stiffness matrix; inverse design; stacking sequence; minimum mass; buckling; Python.')
sub('Highlights')
for s in [
 'Single consistent cylindrical wall solution is used throughout all five numerical studies and the 14-configuration publication comparison.',
 'A 4,128-combination inverse-design screen includes material type, material layering, 3/5/7/9 plies, stacking order and four built-in cellular architectures.',
 'The prescribed full 6×6 stiffness target is enforced by a generalised eigenvalue condition, not individual matrix entries.',
 'Manufacturability, illustrative cost and three-load local-wall buckling are screened separately from minimum stiffness.',
 'Fifteen original computer-generated figures and the complete analysis tables/projects support reproducibility.'
 ]:bullet(s)
fig('F15_workflow.png',1,'Analysis chain from design selection to constrained minimum-mass ranking.')
page()

title('1. Introduction and study scope')
para('Lightweight laminated-composite cellular cores provide multiple interacting design freedoms: periodic wall topology, cell dimensions, constituent elastic constants, fibre orientation, ply count, stacking symmetry and wall thickness. Their combined influence is naturally investigated by energy-equivalent homogenisation, which replaces the periodic wall graph by a continuum with the same macroscopic strain energy. The mathematical and experimental reference used here is Sriharan et al. (2024), Engineering Structures 299, 117106 [1].')
para('The paper-style implementation discretises walls into laminated Timoshenko elements, imposes periodic constraints, solves six independent macroscopic strains and condenses the microstructural strain energy into a symmetric 6×6 stiffness matrix. GUI workspaces and batch calculations share the same Python mechanics backend. This manuscript places emphasis on reproducible engineering studies rather than proposing a new constitutive theory.')
para('Scope decision: ALL numerical studies and recommendations below use wall_element_type="paper" and strip_bending_mode="cylindrical". No free-moment outputs, corrections or cross-formulation rankings are included in any result or design decision. The custom topology example is a separate feasibility demonstration and is not included in the four-family target-C optimisation palette.')
sub('1.1 Research questions')
for x in ['How do stiffness, anisotropy and density scaling vary across built-in periodic geometries?', 'Can a large solved database identify nominal and conservatively screened Pareto-optimal designs?', 'Does the governing wall under compression/shear change with the macroscopic load combination?', 'Can a user-defined graph be mechanically compared with a matched built-in periodic topology?', 'Which combination of built-in geometry, ply count, fibre angles and material arrangement gives the lowest core mass while meeting the complete target C matrix and auxiliary constraints?']:bullet(x)

sub('1.2 Caution about interpretation')
para('The computed stiffnesses use a reduced-order analytical wall model. The plate buckling screen assumes simply supported rectangular walls and classical laminated-plate theory. A lower matrix norm does not imply improved safety, and full homogenised stiffness dominance does not mean every scalar coefficient closely matches the prescribed value. A material-cost index represents relative library weights, not market prices. None of the optimisation results is presented as a mathematically certified global optimum or a substitute for FE eigenbuckling and test validation.')

title('2. Computational formulation and implementation')
sub('2.1 Laminated wall constitutive model')
para('Each wall is defined by an explicit ply stack and evaluated using plane-stress transformed lamina stiffness. The laminate A, B and D matrices follow classical laminate theory, with additional neutral-plane and equivalent transverse-shear operations supplying the cylindrical Timoshenko wall. The published-style model retains the constrained cylindrical wall response for all results reported in this article.')
eq('A = Σ Q̄ₖ (zₖ₊₁ − zₖ)  ;  B = ½Σ Q̄ₖ (z²ₖ₊₁ − z²ₖ)  ;  D = ⅓Σ Q̄ₖ (z³ₖ₊₁ − z³ₖ)')
para('The wall stiffness operator maps local nodal motion and the applied through-thickness strain into element energy. By rotating and assembling the graph elements with physical wall weights, the RVE energy takes the quadratic form shown below. Wall thickness, layup, wall angle and topology affect the assembled coefficients.')
eq('U = ½ uᵀKu + uᵀHε + ½ εᵀRε')
sub('2.2 Periodicity, condensation and engineering properties')
para('Two lattice vectors a₁, a₂ define an in-plane periodic cell. The solver assembles periodic displacement jumps in macroscopic engineering-strain coordinates and removes rigid-body translations with a compatible anchor. The constrained system is solved simultaneously for six right-hand sides. The reduced energy Hessian divided by the periodic RVE volume produces Cᴴ, followed by compliance inversion and scalar elastic-property extraction.')
eq('Cᴴ = V⁻¹ (UᵀKU + UᵀH + HᵀU + R)  ;  Sᴴ = (Cᴴ)⁻¹')
para('In the full internal matrix, the engineering Voigt index order is [11,22,33,23,13,12]. The article-style reduced shell convention instead displays Q44 as 1–3 shear and Q55 as 2–3 shear. Hence the shell output is Q44 = Cᴴ₅₅ and Q55 = Cᴴ₄₄ (1-based conventional indices); no rows or columns are swapped inside the stored full 3D C matrix.')
sub('2.3 Physical quantities')
para('For wall k, let Lₖ, wₖ, tₖ be the element length, physical counting weight and total laminate thickness. At core height h and cell area Acell the relative density, effective mass density and core mass per face area are:')
eq('ρ* = (Σ wₖ Lₖ tₖ) / Acell  ;  ρeff = mcore / (Acell h)  ;  mA = h ρeff')
para('With h=150 mm, a density of 83.088 kg/m³ corresponds to mA=12.463 kg/m². This is core-only mass; the face sheets, adhesive and connections are not included. Material layer density is accumulated by each ply; hybrid walls therefore have mass-weighted material proportions.')
sub('2.4 Software workflow and reproducibility')
para('The implementation uses NumPy/SciPy for matrix assembly and sparse constrained solving, a Python API independent of the PySide6 user interface, JSON input libraries, CSV/Parquet/XLSX batch outputs, and a compact reduced-order buckling evaluator. Reproducible research requires the material and geometry library definitions as well as the case settings; storing the numeric geometry identifiers alone is insufficient.')
fig('F01_builtin_geometries.png',2,'Representative built-in periodic wall graphs and their solver node and wall counts; diagrams are generated from the actual geometry functions, not GUI screenshots.')

page()
title('3. Four built-in geometry and laminate design libraries')
para('The target-C study restricts geometry selection to HEX, SQUARE, MIXED and TRI. Although the previously released Node Explorer supports additional 3–10-node layouts, those exploratory networks are not mixed into the present built-in-only optimisation. Each native graph is parameterised by two geometry variables. Topology is treated as a discrete variable and parameter dimensions are continuously sampled.')
table(['Geometry family','Parameter 1','Parameter 2','Nodes','Walls'],[['HEX','L1 = 35–70 mm','L2 = 35–70 mm','5','3'],['SQUARE','Lx = 35–70 mm','Ly = 35–70 mm','7','4'],['MIXED','Lh = 30–60 mm','R = 1.5–2.5','7','6'],['TRI','Lh = 35–70 mm','R = 0.8–1.2','9','6']], 'Table 1. Built-in geometry ranges used for the new optimisation.')
para('Four material constructions are used: all E-glass/epoxy, all boron/epoxy, an alternating glass–boron–glass family, and an alternating boron–glass–boron family. The hybrid shorthand denotes the outer/centre material ordering pattern and is expanded according to the chosen number of plies; every laminate is symmetric in fibre angles but material sequence is simply alternating and symmetric for odd ply counts. Elastic moduli and densities are taken from the active software library. The cost indices are 1 and 3 per kg for glass/epoxy and boron/epoxy, respectively; these values are illustrative, not supplier quotations.')
table(['Material','E1 (MPa)','E2 (MPa)','G12 (MPa)','ν12','Density (kg/m³)','Relative cost index'],[['E-glass/epoxy','38,600','8,270','4,140','0.26','1,800','1'],['Boron/epoxy','204,000','18,500','5,590','0.23','2,080','3']], 'Table 2. Constituent input values taken from the installed material library.')
para('The ply-count domain is Np ∈ {3, 5, 7, 9}. Angle choices come from {0°, +45°, −45°, 90°}; only selected symmetric sequences are enumerated and their order is retained in the laminate model. Different ply counts therefore change the through-thickness orientation distribution, not simply the label of an otherwise identical homogeneous wall. Each ply is assigned equal thickness within a wall. Any formulation with tply < 0.08 mm is rejected.')
table(['Ply count','Illustrative symmetric arrangements','Explicit material schemes'],[['3','[0/45/0], [90/0/90], [45/−45/45]','Glass, Boron, G/B/G, B/G/B'],['5','[0/90/45/90/0], [0/45/0/45/0]','Same four schemes'],['7','[0/90/45/0/45/90/0]','Same four schemes'],['9','[0/90/45/0/90/0/45/90/0]','Same four schemes']], 'Table 3. Examples from the candidate stacking-sequence generation algorithm; the complete generated palette is provided as machine-readable input.')

page()
title('4. Publication benchmark: 14 laminate configurations')
para('The previously published 14 material/laminate configurations are retained as a source-aligned cylindrical benchmark, with all six reported shell stiffnesses checked independently. This is a numerical consistency test against the published analytical values, not a recalibration of the wall formulation.')
b=bench.copy();b['pct']=100*(b.PaperCyl_MPa-b.Published_MPa)/b.Published_MPa
pr=b.groupby('Property').pct.agg(lambda x:np.max(np.abs(x))).to_dict()
para(f'Across 84 paper reference values, {sum(np.abs(b.pct)<=1)} are within ±1%. The maximum deviation occurs for Q66 ({pr["Q66"]:.3f}%). The other five stiffness terms show much smaller discrepancies. The Q66 difference remains unresolved and is retained as an explicit model-validation limitation when interpreting stiffness-dominated optimisation outcomes.')
table(['Term','Max. absolute publication difference (%)'],[[k,f'{pr[k]:.3f}'] for k in ['Q11','Q22','Q12','Q66','Q44','Q55']], 'Table 4. 14-configuration cylindrical-only paper benchmark errors; signed differences use 100 (Python – paper)/paper.')
fig('F02_journal_benchmark_cyl.png',3,'Cylindrical model differences for 14 published configurations and six reported stiffnesses.')

page()
title('5. Case Study 1 — High-throughput topology–material design atlas')
para('The initial atlas contains 2,400 independent samples distributed equally over HEX, SQUARE, MIXED and TRI. Relative density varies from 0.04 to 0.12, with two pure-material systems and six library stacking sequences. Each design is solved exactly once in the cylindrical wall mode for the present analysis; previously archived alternative-model columns are excluded from every summary and figure.')
med=atlas.groupby('GeometryID').Qmin_paper.median().reindex(['HEX','SQUARE','MIXED','TRI'])
table(['Topology','Designs','Median Qmin (MPa)','Median specific Qmin'],[[g,'600',f'{med[g]:.2f}',f'{atlas[atlas.GeometryID==g].SpecificQmin_paper.median():.3f}'] for g in ['HEX','SQUARE','MIXED','TRI']], 'Table 5. Descriptive cylindrical response statistics for the 2,400-design atlas; independent samples are not pairwise matched across families.')
fig('F03_atlas_specific_stiffness.png',4,'Computed stiffness–density cloud and descriptive specific-stiffness statistics for the cylindrical atlas.')
para('To isolate topology-dependent density scaling, a separate controlled series varies the total wall thickness at fixed geometry/material/layup parameters. A power-law representation of the form Q ~ C (ρ*)ⁿ is fitted to each sequence. Its slope n captures stretching- versus bending-dominated scaling only within the sampled parameter range and model assumptions.')
eq('Qᵢⱼ = Cᵢⱼ (ρ*)ⁿ   ;   n = d log(Qᵢⱼ) / d log(ρ*)')
ff=fits[(fits.Mode=='paper')&(fits.Metric=='Q66_MPa')].groupby('GeometryID').n.median();qq=fits[(fits.Mode=='paper')&(fits.Metric=='Qmin_MPa')].groupby('GeometryID').n.median()
table(['Topology','Median n(Qmin)','Median n(Q66)'],[[g,f'{qq[g]:.3f}',f'{ff[g]:.3f}'] for g in ['HEX','SQUARE','MIXED','TRI']], 'Table 6. Cylindrical-only median power-law exponents across controlled density sweeps.')
fig('F04_density_scaling.png',5,'Fitted density exponents for Qmin and Q66; error bars represent the range across controlled material/layup combinations.')
para('In this model and within the evaluated ranges, Q66 is nearly cubic in relative density for the hexagonal and square families, but nearly linear for the mixed and triangular families. The weaker normal stiffness is approximately linear across all four families. The difference makes topology selection important when shear performance is a governing constraint.')

page()
title('6. Case Study 2 — Robust multi-objective inverse design')
para('The second study reuses the 2,400-case atlas but ranks designs using nominal density, specific normal stiffness, relative material cost and deviation from a target in-plane stiffness ratio Q11/Q22 = 1.5. The earlier design database applies explicit manufacturability rules and a deterministic conservative response scenario rather than claiming probabilistic reliability. All computations and rankings discussed here use cylindrical-only results.')
eq('Minimise [ ρeff, CostIndex, |Q11/Q22 − 1.5|, −Qmin ]')
para('The conservative scenario assumes 5% stiffness variation, 3% density variation and 10% cost variation with a multiplier of 1.65, as defined in the archived four-study experiment. These are hypothetical engineering margins, not a statistically established manufacturing distribution. A candidate is Pareto dominated only when another candidate is no worse in all four criteria and strictly better in at least one.')
table(['Scenario','Cylindrical feasible designs','Non-dominated candidates'],[['Nominal','1,032','79'],['Conservative','932','78']], 'Table 7. Previously computed cylindrical-only multi-objective finite-database screening counts.')
fig('F05_robust_inverse_nominal.png',6,'Nominal stiffness–density and target-ratio distributions used to illustrate inverse-design criteria.')
para('The previously selected trade-off design is database candidate C002152 (TRI), while a separate closest-ratio example achieved approximately Q11/Q22 = 1.497 in the cylindrical model. These results demonstrate inversion against a scalar design target and should not be confused with the full 6×6 target-C optimization developed later in this manuscript.')

page()
title('7. Case Study 3 — Load-dependent wall buckling')
para('For eight previously selected cellular architectures, a seven-load parametric study recovers wall membrane resultants from homogenised macroscopic stresses and evaluates simply supported laminated-wall plate bifurcation. With a single cylindrical interpretation, this represents 8 × 7 = 56 load assessments. The critical-wall identity is determined separately for each macroscopic load vector.')
eq('σ = [σ11, σ22, σ33, τ23, τ13, τ12]ᵀ   ;   λmin = min over walls (λcr,wall)')
fig('F06_buckling_load_interaction.png',7,'Cylindrical-only load-dependent buckling factors and governing-wall counts for eight selected architectures.')
para('In the archived MIXED example C001592, the uniaxial compression multiplier is approximately 8.415 (wall 6), whereas the specified compressive/shear interaction reduces the multiplier to about 0.647 (wall 4). The governing wall changes as the load direction changes. This screening result explains why a stiffness-optimal core cannot be accepted without considering the intended loading.')
para('The buckling implementation uses free-membrane ABD condensation when required and a simply supported multi-mode Ritz/Navier plate basis. Mode pairs reported below are dominant terms of that approximate basis, not full-field finite-element eigenmodes. The core geometry, wall junctions, facesheet support and imperfections are not represented in a three-dimensional stability analysis.')

page()
title('8. Case Study 4 — Custom periodic topology and laminate tailoring')
para('A custom six-wall diagonal-braced periodic plug-in was constructed using node coordinates, connectivity, lattice vectors and wall weights. It was evaluated against a matched built-in square control for 120 design settings. Matching ensures comparable geometry parameters, constituent material, layup and relative density; only the graph connections differ. In this section, only the 120 cylindrical pairs (240 solves) are included.')
cm=custom[custom.Mode=='paper'];delta=cm.Qmin_MPa_Braced_vs_Square_pct
para(f'Across all 120 matched cylindrical comparisons, the diagonal-braced layout lowers Qmin relative to its square control; the median change is {delta.median():+.2f}%. In contrast, Q66 increases substantially, illustrating the potential to tailor anisotropy and shear resistance by changing topology. The custom plug-in remains outside the four-family target-C selection palette.')
fig('F07_custom_matched_cyl.png',8,'Matched weak-direction normal and in-plane shear stiffness of custom diagonal-braced and square reference networks.')
para('The supplied plug-in is a structural graph definition that can be registered in the Geometry Library. Its successful Python-side discovery and analysis were tested. A screenshot of an operational GUI editing session is intentionally not claimed, because the PySide6 window was not launched during this calculation.')

page()
title('9. Case Study 5 — Minimum-weight design for a prescribed full stiffness matrix')
sub('9.1 Target and physically meaningful constraint')
para('The full target stiffness is fixed to 80% of the solver’s published configuration-2 cylindrical benchmark. The matrix is symmetric and positive definite and refers to a three-dimensional orthotropic continuum in MPa with internal engineering Voigt order [11,22,33,23,13,12]. No shell-output index swap is applied when assessing matrix adequacy. The target is fixed before running the optimisation.')
table(['','11','22','33','23','13','12'],[[['11','22','33','23','13','12'][i]]+[f'{T[i,j]:.3f}' for j in range(6)] for i in range(6)],'Table 8. Prescribed target C matrix, in MPa. All rows/columns are in the solver’s internal 3D engineering Voigt order.')
para('The objective is to minimise core-only areal mass at h=150 mm. Rather than demanding exact coefficient-by-coefficient equality, target satisfaction requires the achieved stiffness tensor to dominate the target for every possible macroscopic engineering strain. This is tested by the minimum generalised eigenvalue of the achieved and target 6×6 matrices. A 5% margin on the target is required in the new study.')
eq('minimise mA = h ρeff     such that     λmin(Cᴴ, Ctarget) ≥ 1.05')
para('A value λmin ≥ 1.05 is equivalent to Cᴴ − 1.05 Ctarget being positive semidefinite in the symmetric numerical setting. It does NOT require the Frobenius matrix mismatch ||Cᴴ − Ctarget||F / ||Ctarget||F to be small. Both quantities are therefore reported independently.')
sub('9.2 Design variables and fixed bounds')
table(['Class','Decision or constraint','Search definition'],[['Discrete geometry','Built-in topology','HEX, SQUARE, MIXED, TRI'],['Geometry','Two continuous parameters','Bounds in Table 1, six seeded sets per family'],['Material and ply material','4 constituent patterns','Glass, Boron, G/B/G, B/G/B'],['Layup','Number of plies','3, 5, 7 or 9'],['Fibre angles','Symmetric discrete sequences','0°, +45°, −45°, 90°'],['Total wall thickness','Continuous density scaling','0 < ρ* ≤ 0.15'],['Manufacturing','Minimum wall and individual ply thickness','twall ≥ 0.35 mm; tply ≥ 0.08 mm'],['Manufacturing','Largest/smallest wall thickness','≤ 1.9'],['Illustrative material budget','Index per m²','≤ 50 relative cost units'],['Stiffness','Generalised eigenvalue','≥ 1.05'],['Stability','Uniaxial, biaxial, combined stress load cases','Each λbuckling ≥ 1.5']], 'Table 9. Variable domain and separate engineering constraints in the cylindrical-only inverse design.')
para('The load cases use macroscopic stress vectors in MPa: (i) [−1, 0, 0, 0, 0, 0], (ii) [−1, −0.5, 0, 0, 0, 0] and (iii) [−1, −0.5, 0, 0, 0, +0.25]. The cost ceiling is a deliberately illustrative constraint for comparing materials, not a purchase or production cost limit. The fixed geometry is prismatic with no face sheets or adhesive.')
sub('9.3 Candidate-generation and optimisation protocol')
para('The seeded grid generates 4,128 candidate records. Every record uses a built-in geometry function, an explicit laminate stack and the same paper-style wall element. First, candidates are screened at relative density 0.15. Promising profiles then undergo uniform wall-thickness reduction using binary searches for the 1.05 stiffness threshold subject to ply and wall thickness. A material-balanced and topology-balanced shortlist is further refined under the three buckling loads and illustrative cost ceiling. The buckling basis uses up to five modes per direction during ranking; selected saved designs are recalculated with six modes per direction.')
para('The selected sampling, discrete stacking palette and follow-up optimisation are finite. In the reported workflow the 1,177 stiffness refinements include 958 prioritised records plus 219 auxiliary records from an earlier exploratory pass. A shortlist of 136 is screened for buckling, of which 110 pass all declared checks. This is a computational finite-pool optimum, not proof of a global minimum over all allowed layer arrangements or continuous geometry.')
fig('F08_search_funnel.png',9,'The finite search pipeline and original full-C screening response at maximum allowed relative density.')

page()
sub('9.4 Optimal geometry, material and ply count')
para(f'The lightest passing shortlisted solution is {best.ID}: {best.Geometry} topology, {best.MaterialScheme} material, {int(best.PlyCount)} plies [0/45/0], with L1 = L2 = 52.5 mm. Uniform wall thickness refinement gives average density {best.Density:.4f} kg/m³, areal mass {best.ArealMass:.4f} kg/m² and relative density {best.RelDensity:.6f}. The minimum wall thickness is {best.MinThickness:.4f} mm and the corresponding minimum ply thickness is {best.MinPlyThickness:.4f} mm.')
fig('F09_geometry_material_ply_choice.png',10,'Minimum passing density by built-in geometry and heatmap of material arrangement versus ply count (only shortlisted buckling-screened profiles).')
fig('F12_selected_geometry_layup.png',11,'Selected paper-style cylindrical hexagonal graph and three-ply [0/45/0] wall stacking sequence.')
minby=ok.groupby('Geometry').Density.min().reindex(['HEX','SQUARE','MIXED','TRI'])
table(['Geometry','Lightest screened feasible density (kg/m³)','Areal mass (kg/m²)'],[[g,f'{minby[g]:.2f}',f'{minby[g]*.15:.2f}'] for g in ['HEX','SQUARE','MIXED','TRI']], 'Table 10. Minimum density by geometry in the 136-profile screened shortlist; no global-rank claim is implied.')
sub('9.5 A lower-cost near-weight alternative')
alt=ok[ok.ID=='P00148'].iloc[0]
para(f'The five-ply symmetric hybrid hexagonal candidate {alt.ID}, with alternating boron/glass layers and [0/90/45/90/0] fibre orientations, weighs {alt.Density:.3f} kg/m³ compared with {best.Density:.3f} kg/m³ for the pure-boron minimum. Its core areal mass is {alt.ArealMass:.3f} kg/m² and illustrative cost index is {alt.TotalCostPerM2:.3f}/m² versus {best.TotalCostPerM2:.3f}/m² for the minimum. This gives nearly the same mass with a materially lower cost index, subject to the assumed constituent costs and the current buckling model.')
table(['ID','Geometry','Material','Plies / angle sequence','Density kg/m³','Mass kg/m²','Cost index/m²'],[[r.ID,r.Geometry,r.MaterialScheme,str(int(r.PlyCount))+' / '+str(r.Angles).replace(', ','/').replace('[','[').replace('] ',']'),f'{r.Density:.2f}',f'{r.ArealMass:.2f}',f'{r.TotalCostPerM2:.2f}'] for _,r in ok.nsmallest(6,'Density').iterrows()], 'Table 11. Six lowest-mass passing designs in the screened shortlist.')

page()
sub('9.6 Achieved stiffness matrices and fidelity')
para(f'The minimum-mass solution has λmin(Cᴴ,Ctarget) = {best.Adequacy:.6f}. Thus the complete stiffness tensor clears the required 1.05 lower-bound margin. However, the relative matrix Frobenius gap is {best.MatrixGap*100:.1f}%; several coefficients exceed their targets substantially. The hybrid alternative reduces that mismatch to {alt.MatrixGap*100:.1f}% while remaining within the stiffness and buckling constraints. The distinction between energetic adequacy and coefficient matching is decisive.')
table(['','11','22','33','23','13','12'],[[['11','22','33','23','13','12'][i]]+[f'{C[i,j]:.2f}' for j in range(6)] for i in range(6)],'Table 12. Achieved full Cᴴ of design P00050 in MPa, cylindrical-only solution, internal engineering Voigt order.')
fig('F11_target_and_achieved_matrices.png',12,'Target tensor and the achieved tensors for the minimum-mass and lower-cost hybrid solutions. Consistent colour scale (MPa).')
fig('F10_mass_cost_matrix_matching.png',13,'Full-matrix discrepancy versus density and relative material cost versus areal mass; the minimum satisfying design need not be the closest coefficient-by-coefficient match.')
para('For applications that require direct tensor matching, a second-stage objective should minimise the full-matrix difference, possibly with user-specified tolerances on individual coefficients and off-diagonal couplings. Such a matching objective can conflict with minimising density and is separate from the present lower-bound target-stiffness design problem.')

page()
sub('9.7 Local wall stability, governing walls and modes')
para(f'For design {best.ID}, the six-mode verification gives buckling factors {best.Buck_uniaxial:.3f}, {best.Buck_biaxial:.3f} and {best.Buck_combined:.3f} under the prescribed uniaxial, biaxial and combined loads, respectively. The critical wall is wall {int(best.Wall_uniaxial)} for uniaxial loading, wall {int(best.Wall_biaxial)} for biaxial loading and wall {int(best.Wall_combined)} for combined loading. The dominant Ritz term is (1,1) for each load. The wall-buckling factors are preliminary and include no imperfections or FE-validated junction restraint.')
fig('F13_buckling_selected.png',14,'Buckling multipliers, required factor 1.5 and governing-wall changes for four representative feasible designs.')
fig('F14_constraint_margins.png',15,'Ratio of achieved to required constraint values for design P00050. Ratios ≥ 1 indicate the stated screening rule is met.')
table(['Verification','Requirement','P00050 obtained','Result'],[['Full-C adequacy','≥1.05',f'{best.Adequacy:.4f}','Pass'],['Effective density','ρ*≤0.15',f'{best.RelDensity:.5f}','Pass'],['Wall thickness','≥0.35 mm',f'{best.MinThickness:.3f} mm','Pass'],['Individual ply thickness','≥0.08 mm',f'{best.MinPlyThickness:.3f} mm','Pass'],['Wall thickness ratio','≤1.90',f'{best.ThicknessRatio:.3f}','Pass'],['Illustrative cost','≤50 /m²',f'{best.TotalCostPerM2:.3f}','Pass'],['Uniaxial buckling','≥1.50',f'{best.Buck_uniaxial:.3f}','Pass'],['Biaxial buckling','≥1.50',f'{best.Buck_biaxial:.3f}','Pass'],['Combined buckling','≥1.50',f'{best.Buck_combined:.3f}','Pass']], 'Table 13. All declared cylindrical-only constraints for the selected minimum-mass design.')

page()
title('10. Discussion, verification and limitations')
sub('10.1 What the five studies establish')
para('The first atlas demonstrates that geometry and laminate choices affect stiffness scaling and specific mechanical response. The multi-objective study demonstrates that a single nominal stiffness optimum is not necessarily a strong inverse-design solution once density, cost, ratio and manufacturing margins enter. The buckling interaction study establishes that different load vectors can activate different wall weaknesses. The custom graph confirms that topology can be supplied through a reusable API. The new inverse-design study integrates discrete and continuous choices with a full six-dimensional stiffness requirement and auxiliary safety screens.')
sub('10.2 Cylinder-only scope and validation boundary')
para('Every stiffness calculation in the primary studies and new search uses the original paper-style cylindrical wall solver. Results from anisotropic free-moment options are deliberately omitted from the final manuscript. The 14-configuration article comparison is retained as an independent benchmark and its unresolved Q66 difference (maximum ~14.54%) must be addressed before quantitative claims of full model validation can be extended to new topologies and laminate stacks. The benchmark does not constitute FE validation of the new selected optimum.')
sub('10.3 Buckling, manufacturing and material limitations')
para('The local-wall buckling constraint uses a linear simply supported laminated-plate model, not the stability of a joined cellular sandwich panel. It does not explicitly account for imperfect geometry, adhesive joints, residual stresses, material failure, damage progression, skin restraint, short-wave local modes beyond the chosen Ritz basis, or global sandwich-panel buckling. The ply thickness and wall thickness limits are illustrative fabrication parameters and must be replaced with process-specific limits for a real product. The cost ceiling is an assumed relative index; material prices and processing costs have not been verified.')
sub('10.4 Finite search and design-fidelity limitation')
para('The minimum found is the lowest density in a screen of 136 selected profiles after an initial 4,128-combination exploration and selected thickness refinements. It is not a proof of optimality in the continuous design space. Uniform thickness scaling is applied in the search, although separate wall-thickness factors can be introduced in a later constrained optimisation. The selected design satisfies Cᴴ≽1.05Ctarget but does not reproduce the target entries; a near-target objective would yield a different Pareto frontier. Both sources of approximation should appear whenever the design is cited.')
sub('10.5 Reproducibility and recommended verification')
para('The research archive contains the target matrix, 4,128 candidate definitions, the stiffness and buckling result tables, selected full stiffness matrices and loadable .rvejson project files. It includes a verification utility that reloads saved projects, solves them using the same cylindrical model and compares their 6×6 matrices and screened stability margins. The companion version of the software is not required to load the figures, but is required to recalculate physical answers.')
para('A defensible next validation phase should independently generate periodic 3D shell or solid FE models for P00050 and P00148, confirm the effective C tensor with compatible PBCs, and undertake wall/panel eigenbuckling with realistic joint boundary conditions. Physical fabrication tolerances and actual cost models should then replace illustrative constraints before engineering deployment.')

title('11. Conclusions')
for x in [
 'A consistent paper-style cylindrical-bending mechanics pathway can support graph-based homogenisation, batch exploration, inverse design and local-wall load recovery without mixing incompatible wall models.',
 'A controlled 2,400-design study reveals substantially different Q66 density-scaling exponents across the four built-in geometry families, while the normal-stiffness scaling is close to linear.',
 'A new 4,128-combination search over geometry, two constituent systems, hybrid material arrangements, four ply counts and symmetric stacking sequences produced 110 passing candidates from the 136-profile buckling shortlist.',
 f'The lightest passing candidate is a [0/45/0] three-ply boron/epoxy HEX core with ρeff={best.Density:.2f} kg/m³, mA={best.ArealMass:.2f} kg/m², generalised stiffness factor {best.Adequacy:.4f} and combined-load buckling multiplier {best.Buck_combined:.3f}.',
 f'A five-ply hybrid HEX alternative provides comparable density ({alt.Density:.2f} kg/m³), with lower assumed cost ({alt.TotalCostPerM2:.2f} rather than {best.TotalCostPerM2:.2f} relative cost units/m²).',
 'The study establishes a minimum required-stiffness solution, not an entrywise match of the target C matrix. The Q66 publication mismatch, finite-candidate nature and reduced-order buckling model remain explicit limits.'
 ]:bullet(x)

page()
title('Appendix A. Reproduction protocol and input files')
para('Install the supplied RVE Studio v1 package with NumPy, SciPy, pandas and the optional plotting stack. The reproduction package includes a data-specific cylindrical search script, JSON records of geometry/material/ply candidates and the exact numerical stiffness target. All source data used to regenerate article figures are exported with paper-style cylindrical columns only. The saved RVE projects can be reopened with the standard API or Single RVE Designer.')
para('Recommended commands (after installing RVE Studio and extracting the research ZIP):')
p=DOC.add_paragraph();r=p.add_run('python verify_cylindrical_results.py\npython run_cylindrical_study.py screen\n# Full refinement is performed in documented numbered batches.\npython make_figures.py');r.font.name='Consolas';r.font.size=Pt(8.5)
para('For the main minimum-mass design, the saved project is results/project_P00050.rvejson and the achieved matrix is results/C_achieved_P00050_MPa.csv. All dimensions are in mm and stiffnesses are in MPa. The fixed target comes from the source published configuration 2 through the cylindrical solver; it is supplied as results/C_target_MPa.csv to avoid ambiguity.')
sub('Appendix A.1 Audit of the 14 published configurations')
prop=['Q11','Q22','Q12','Q66','Q44','Q55']
rows=[]
for n in sorted(bench.Config.unique()):
 d=bench[bench.Config==n].set_index('Property');rows.append([n]+[f'{100*(d.loc[k,"PaperCyl_MPa"]/d.loc[k,"Published_MPa"]-1):+.2f}' for k in prop])
table(['Cfg']+prop,rows,'Table A1. Signed percentage differences for each of 14 configurations versus published paper analytical values, cylindrical-only implementation. Values are percent.')
sub('Appendix A.2 Key output files')
for name in ['initial_screen.csv: initial 4,128 designed laminate/material/geometry configurations at ρ*=0.15.','stiffness_refined.csv: 1,177 continuous wall-thickness refinement records.','buckling_refined.csv: 136 buckling-screened candidates, including 110 that pass all declared checks.','C_target_MPa.csv: six-by-six prescribed target in internal engineering Voigt ordering.','C_achieved_P00050_MPa.csv and project_P00050.rvejson: selected optimum output and loadable project.','C_achieved_P00148_MPa.csv: lower-cost hybrid alternative.','paper_only_*.csv: prior four-study outputs filtered to cylindrical results only.','figures/F01–F15: original figure images embedded in this report.']:bullet(name)

title('References')
para('[1] Sriharan, J., Dias, M., Fernando, D. and Adhikari, S. (2024). “A unified approach for the prediction of the effective properties of laminated composite cellular core.” Engineering Structures, 299, 117106. https://doi.org/10.1016/j.engstruct.2023.117106.')
para('[2] RVE Studio v1.0.0 software release. Packaged source files rvestudio/models.py, laminate.py, walls.py, geometry.py, homogenization.py, properties.py, physical.py, load_recovery.py and buckling.py. Accompanying reproducibility archive, 30 October 2026.')
para('[3] RVE Studio v1 previous four-case-study numerical dataset. Original seeded 2,400-design atlas, robust screening and matched custom-geometry results, archived alongside the current case-study computations. Cylindrical columns only are used in this manuscript.')
DOC.core_properties.title='RVE Studio v1 — Cylindrical-Bending Five Case Studies and Target C Inverse Design'
DOC.core_properties.subject='Cylindrical-only RVE homogenisation; constrained minimum mass and ply selection'
DOC.core_properties.keywords='homogenisation, cylindrical bending, target C matrix, plies, inverse optimisation'
DOC.save(OUT)
print('WROTE',OUT,OUT.stat().st_size,'figures',15,'best',best.ID,best.Density)
