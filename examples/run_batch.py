from rvestudio.database import run_batch

run_batch('RVE_5000_python.csv','RVE_5000_python_results.csv',workers=1,case_limit=20,
          progress=lambda i,n,r: print(f"[{i}/{n}] {r['CaseID']} {r['Status']}"))
