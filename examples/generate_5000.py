from rvestudio.database import generate_parametric_database

df=generate_parametric_database('RVE_5000_python.csv',num_cases=5000,seed=20260906)
print(df.GeometryID.value_counts())
print('Generated',len(df),'cases')
