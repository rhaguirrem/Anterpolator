import os
import struct
import glob

def analyze_bmf(filepath):
    try:
        filesize = os.path.getsize(filepath)
        basename = os.path.basename(filepath)
        
        with open(filepath, 'rb') as f:
            data = f.read(131072)
        
        u32_offsets = range(8, 128, 4)
        u64_offsets = [40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120]
        
        u32_values = {}
        for off in u32_offsets:
            if off + 4 <= len(data):
                u32_values[off] = struct.unpack('<I', data[off:off+4])[0]
                
        u64_values = {}
        for off in u64_offsets:
            if off + 8 <= len(data):
                u64_values[off] = struct.unpack('<Q', data[off:off+8])[0]
                
        strings = ['created', 'description', 'lower_x', 'modified']
        str_offsets = {}
        for s in strings:
            idx = data.find(s.encode('ascii'))
            str_offsets[s] = idx if idx != -1 else "Not Found"
            
        return {
            'basename': basename,
            'size': filesize,
            'u32': u32_values,
            'u64': u64_values,
            'strings': str_offsets
        }
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

folder = r'C:\Users\raguirre\OneDrive - SRK Consulting\Documents\Trabajo\SRK\Misc\BMF Examples'
files = glob.glob(os.path.join(folder, '*.bmf'))

if not files:
    print(f"No .bmf files found in {folder}")
else:
    results = []
    for f in files:
        res = analyze_bmf(f)
        if res:
            results.append(res)

    for r in results:
        print(f"File: {r['basename']} | Size: {r['size']}")
        u32_str = ", ".join([f"{off}:{val}" for off, val in sorted(r['u32'].items())])
        print(f"  u32: {u32_str}")
        u64_str = ", ".join([f"{off}:{val}" for off, val in sorted(r['u64'].items())])
        print(f"  u64: {u64_str}")
        print(f"  Strings: {r['strings']}")
        print("-" * 20)

    if results:
        u32_keys = sorted(results[0]['u32'].keys())
        u64_keys = sorted(results[0]['u64'].keys())
        
        print("\nSummary Analysis:")
        constant_u32 = [k for k in u32_keys if all(r['u32'].get(k) == results[0]['u32'].get(k) for r in results)]
        print(f"Constant u32 offsets: {constant_u32}")
        
        size_matches_u32 = [k for k in u32_keys if all(r['u32'].get(k) == r['size'] for r in results)]
        print(f"u32 offsets matching filesize: {size_matches_u32}")

        size_matches_u64 = [k for k in u64_keys if all(r['u64'].get(k) == r['size'] for r in results)]
        print(f"u64 offsets matching filesize: {size_matches_u64}")

        for off in [48, 112, 120]:
            match_count = sum(1 for r in results if r['u64'].get(off) == r['size'])
            vals = [r['u64'].get(off) for r in results]
            print(f"Offset {off} (u64) matches filesize in {match_count}/{len(results)} files. Values: {vals}")
