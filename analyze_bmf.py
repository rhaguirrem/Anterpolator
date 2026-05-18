import os, struct, glob
folder = r'C:\Users\raguirre\OneDrive - SRK Consulting\Documents\Trabajo\SRK\Misc\BMF Examples'
files = glob.glob(os.path.join(folder, '*.bmf'))
print(f'Found {len(files)} files.')
results = []
for fpath in files:
    size = os.path.getsize(fpath)
    with open(fpath, 'rb') as f: data = f.read(131072)
    u32s = {o: struct.unpack('<I', data[o:o+4])[0] for o in range(8, 128, 4) if o+4 <= len(data)}
    u64s = {o: struct.unpack('<Q', data[o:o+8])[0] for o in [40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120] if o+8 <= len(data)}
    strs = {s: data.find(s.encode('ascii')) for s in ['created','description','lower_x','modified']}
    print(f'File: {os.path.basename(fpath)} | Size: {size}')
    print(f'  u32: {u32s}')
    print(f'  u64: {u64s}')
    print(f'  Strings: {strs}')
    results.append({'size': size, 'u32': u32s, 'u64': u64s})
if results:
    print('\nSummary Analysis:')
    for off in range(8, 128, 4):
        vals = [r['u32'].get(off) for r in results]
        if all(v == vals[0] for v in vals): print(f'u32 offset {off} is constant: {vals[0]}')
        if all(v == r['size'] for v, r in zip(vals, results)): print(f'u32 offset {off} matches filesize')
    for off in [40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120]:
        vals = [r['u64'].get(off) for r in results]
        if all(v == vals[0] for v in vals): print(f'u64 offset {off} is constant: {vals[0]}')
        if all(v == r['size'] for v, r in zip(vals, results)): print(f'u64 offset {off} matches filesize')
