import bmf_standalone_exporter as bmf
bmf_file = "test_output.bmf"
table = bmf.load_bmf_table(bmf_file)
print(f"Table keys: {list(table.keys())}")
if "data" in table:
    print(f"Data type: {type(table['data'])}")
    if hasattr(table['data'], 'shape'):
        print(f"Data shape: {table['data'].shape}")
