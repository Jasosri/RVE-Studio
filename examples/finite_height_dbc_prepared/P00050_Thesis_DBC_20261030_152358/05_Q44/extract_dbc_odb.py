from __future__ import print_function
import json, os, sys, traceback
from odbAccess import openOdb
odb=None
out_path='case_results.json'
odb_path=sys.argv[-2] if len(sys.argv)>=3 else sys.argv[-1]
if len(sys.argv)>=3: out_path=sys.argv[-1]
def repo_get(repo,name):
    try:return repo[name]
    except Exception:pass
    try:
        target=str(name).upper()
        for key in repo.keys():
            if str(key).upper()==target:return repo[key]
    except Exception:pass
    return None
def write(data):
    with open(out_path,'w') as f: json.dump(data,f,indent=2)
def find_allse(step):
    vals=[]
    for rk in step.historyRegions.keys():
        reg=step.historyRegions[rk]
        for hk in reg.historyOutputs.keys():
            if str(hk).upper()=='ALLSE':
                data=reg.historyOutputs[hk].data
                if data:vals.append((str(rk),float(data[-1][1])))
    if not vals: raise RuntimeError('ALLSE history output not found')
    vals.sort(key=lambda x:abs(x[1]),reverse=True)
    return vals[0][1],vals
try:
    odb=openOdb(path=odb_path,readOnly=True)
    step=repo_get(odb.steps,'DBC')
    if step is None: raise RuntimeError('DBC step not found')
    energy,candidates=find_allse(step)
    result={'status':'COMPLETED','odb':odb_path,'step':'DBC','strain_energy_ALLSE':float(energy),'allse_candidates':candidates,'frames':len(step.frames)}
    write(result); print(json.dumps(result,indent=2))
except Exception as exc:
    err={'status':'ERROR','odb':locals().get('odb_path',''),'error':str(exc),'traceback':traceback.format_exc()}
    try:write(err)
    except Exception:pass
    print(json.dumps(err,indent=2)); sys.exit(2)
finally:
    if odb is not None:
        try:odb.close()
        except Exception:pass
