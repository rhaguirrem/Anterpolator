import os
import struct

base_dir = r"C:\Users\raguirre\OneDrive - SRK Consulting\Documents\Trabajo\SRK\Misc\BMF Examples"
targets = ["created", "description", "lower_x", "modified"]

results = []
for f in os.listdir(base_dir):
    if f.lower().endswith(".bmf"):
        path = os.path.join(base_dir, f)
        size = os.path.getsize(path)
        with open(path, "rb") as bmf:
            header = bmf.read(128)
            bmf.seek(0)
            scan_data = bmf.read(128 * 1024)
            
            u32 = [struct.unpack("<I", header[i:i+4])[0] for i in range(0, 128, 4)]
            u64 = [struct.unpack("<Q", header[i:i+8])[0] for i in range(0, 128, 8)]
            
            found = {t: scan_data.find(t.encode()) for t in targets}
            results.append({"name": f, "size": size, "u32": u32, "u64": u64, "found": found})

for r in results:
    print(f"File: {r['name']} Size: {r['size']}")
    print(f"  Watch: 32:{r['u32'][8]}/{r['u64'][4]}, 40:{r['u32'][10]}/{r['u64'][5]}, 48:{r['u32'][12]}/{r['u64'][6]}, 56:{r['u32'][14]}/{r['u64'][7]}, 64:{r['u32'][16]}/{r['u64'][8]}, 72:{r['u32'][18]}/{r['u64'][9]}, 80:{r['u32'][20]}/{r['u64'][10]}, 88:{r['u32'][22]}/{r['u64'][11]}, 96:{r['u32'][24]}/{r['u64'][12]}, 104:{r['u32'][26]}/{r['u64'][13]}, 112:{r['u32'][28]}/{r['u64'][14]}, 120:{r['u32'][30]}/{r['u64'][15]}")
    print(f"  Strings: {r['found']}")

if results:
    print("\nSUMMARY:")
    for off in range(0, 128, 4):
        v32s = [r["u32"][off//4] for r in results]
        if all(v == r["size"] for v, r in zip(v32s, results)): print(f"SIZE CORRELATION: u32@{off}")
        if all(v == v32s[0] for v in v32s): print(f"CONSTANT: u32@{off}={v32s[0]}")
    for off in range(0, 128, 8):
        v64s = [r["u64"][off//8] for r in results]
        if all(v == r["size"] for v, r in zip(v64s, results)): print(f"SIZE CORRELATION: u64@{off}")
