import csv
import anterpolator3DViewer as viewer
import bmf_standalone_exporter as bmf

csv_file = "test_input.csv"
bmf_file = "test_output.bmf"
data = [
    ["# Comment line"],
    ["x", "y", "z", "grade", "domain"],
    ["1.0", "2.0", "3.0", "0.5", "10"],
    ["4.0", "5.0", "6.0", "0.8", "20"]
]

with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter=";")
    writer.writerows(data)

viewer.export_csv_grid_to_bmf(
    csv_file, 
    bmf_file, 
    delimiter=";", 
    header_line=2, 
    x_col="x", 
    y_col="y", 
    z_col="z", 
    value_cols=["grade", "domain"]
)

table = bmf.load_bmf_table(bmf_file)
df = table.get("dataframe")

# Report round-trip
expected_count = 2
expected_cols = ["x", "y", "z", "grade", "domain"]

if df is not None:
    row_count = len(df)
    cols = list(df.columns)
    success = row_count == expected_count and all(c in cols for c in expected_cols)
    print(f"Rows count: {row_count}")
    print(f"Columns: {cols}")
    print(f"Success: {success}")
else:
    print("Success: False (No dataframe found)")
