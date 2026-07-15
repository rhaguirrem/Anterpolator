import csv
import anterpolator3DViewer as viewer
import bmf_standalone_exporter as bmf
import pandas as pd

csv_file = "test_input.csv"
bmf_file = "test_output.bmf"
data = [
    ["# Comment line"],
    ["x", "y", "z", "grade", "domain"],
    ["1.0", "11.0", "21.0", "0.5", "10"],
    ["2.0", "12.0", "22.0", "0.8", "20"]
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
    # Check if the values actually match the input data
    # (Filtering out any grid-based padding if present)
    actual_data = df[df['grade'] > 0]
    row_count = len(actual_data)
    cols = list(df.columns)
    
    # Check specifically for the values we put in
    val_check = (1.0 in actual_data['x'].values) and (0.8 in actual_data['grade'].values)
    
    print(f"Total rows in BMF: {len(df)}")
    print(f"Data rows (grade > 0): {row_count}")
    print(f"Columns: {cols}")
    print(f"Value check passed: {val_check}")
else:
    print("No dataframe found")
