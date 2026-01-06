import pandas as pd

# SCL
scl = pd.read_csv('scl_results_fine.csv')
scl['label'] = scl['Composition'].astype(str) + ' (' + scl['Source'].astype(str) + ')'
scl_list = sorted(scl['label'].dropna().unique().tolist())
with open('scl_options.txt','w',encoding='utf-8') as f:
    f.write('\n'.join(scl_list))

# MSTDB
mstdb = pd.read_csv('Molten_Salt_Thermophysical_Properties.csv')
labels = []
for _, row in mstdb.iterrows():
    formula = str(row.get('Formula','')).strip()
    comp = str(row.get('Composition (Mole %)','')).strip()
    if not formula:
        continue
    if '-' not in comp or comp == 'Pure Salt':
        labels.append(formula)  # pure salt
    else:
        try:
            percents = [p.strip() for p in comp.split('-')]
            parts = [p.strip() for p in formula.split('-')]
            if len(percents) == len(parts):
                labels.append('-'.join(f'{perc}{part}' for perc, part in zip(percents, parts)))
        except Exception:
            pass
mstdb_list = sorted(set(labels))
with open('mstdb_options.txt','w',encoding='utf-8') as f:
    f.write('\n'.join(mstdb_list))

# Experimental sources
meas = pd.read_excel('TC_Measurement_Data.xlsx')
srcs = sorted(meas['Source'].dropna().unique().tolist())
with open('experimental_sources.txt','w',encoding='utf-8') as f:
    f.write('\n'.join(srcs))

print('Wrote: scl_options.txt, mstdb_options.txt, experimental_sources.txt')