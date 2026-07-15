import csv
import anterpolator3DViewer as viewer
import bmf_standalone_exporter as bmf

csv_file = "test_input.csv"
bmf_file = "test_output.bmf"
data = [
    ["# Comment line"],
    ["x", "y", "z", "grade", "domain"],
    ["1.0", "11.0", "21.0", "0.5", "10"],
    ["2.0", "12.0", "23.0", "0.8", "20"]
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

if df is not None:
    # Filter for non-zero grades to find our input points
    # Since it is a grid export, original points are preserved in the grid
    data_rows = df[df['grade'] > 0]
    row_count = len(data_rows)
    cols = list(df.columns)
    
    # Verify input values are present
    has_p1 = ((df['x'] == 1.0) & (df['y'] == 11.0) & (df['grade'] == 0.5)).any()
    has_p2 = ((df['x'] == 2.0) & (df['y'] == 12.0) & (df['grade'] == 0.8)).any()
    
    success = row_count == 2 and has_p1 and has_p2 and all(c in cols for c in ["x", "y", "z", "grade", "domain"])
    
    print(f"Data Rows: {row_count}")
    print(f"Columns: {cols}")
    print(f"Values Match: {has_p1 and has_p2}")
    print(f"Success: {success}")
else:
    print("Success: False (No dataframe found)")
