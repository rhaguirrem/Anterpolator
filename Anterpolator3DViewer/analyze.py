import os
import struct

base_dir = r"C:\Users\raguirre\OneDrive - SRK Consulting\Documents\Trabajo\SRK\Misc\BMF Examples"
targets = ["created", "description", "lower_x", "modified"]

if not os.path.exists(base_dir):
    print(f"Error: {base_dir} not found")
    exit()

files = [f for f in os.listdir(base_dir) if f.lower().endswith(".bmf")]
if not files:
    print("No .bmf files found.")
    exit()

results = []
for f in files:
    path = os.path.join(base_dir, f)
    size = os.path.getsize(path)
    with open(path, "rb") as bmf:
        header = bmf.read(128)
        bmf.seek(0)
        scan_data = bmf.read(min(size, 128 * 1024))
        
        u32 = [struct.unpack("<I", header[i:i+4])[0] for i in range(0, 128, 4)]
        u64 = [struct.unpack("<Q", header[i:i+8])[0] for i in range(0, 128, 8)]
        
        found = {t: scan_data.find(t.encode()) for t in targets}
        results.append({"name": f, "size": size, "u32": u32, "u64": u64, "found": found})

for r in results:
    print(f"File: {r['name']} Size: {r['size']}")
    print(f"  Strings: {r['found']}")

print("\nSUMMARY:")
for off in range(0, 128, 4):
    v32s = [r["u32"][off//4] for r in results]
    if all(v == r["size"] for v, r in zip(v32s, results)): print(f"SIZE CORRELATION: u32@{off}")
    if all(v == v32s[0] for v in v32s): print(f"CONSTANT: u32@{off}={v32s[0]}")
for off in range(0, 128, 8):
    v64s = [r["u64"][off//8] for r in results]
    if all(v == r["size"] for v, r in zip(v64s, results)): print(f"SIZE CORRELATION: u64@{off}")
