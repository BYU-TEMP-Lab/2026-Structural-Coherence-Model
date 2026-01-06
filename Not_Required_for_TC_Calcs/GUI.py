import tkinter as tk
from tkinter import IntVar
import pandas as pd

def show_selected():
    selected_indices = [i for i, var in enumerate(checkbox_vars) if var.get()]
    selected_items = [compounds[i] for i in selected_indices]

    # Clear the frame before updating with new selections
    for widget in display_frame.winfo_children():
        widget.destroy()

    # Create Entry widgets only for the selected items
    entry_widgets = []
    for i, item in enumerate(selected_items):
        entry = tk.Entry(display_frame, width=10)
        label = tk.Label(display_frame, text=f"{item}")
        entry.grid(row=i, column=0)
        label.grid(row=i, column=1)
        entry_widgets.append(entry)

    # Store the entry widgets in the selected_compounds list
    selected_compounds.extend(entry_widgets)

def save_data():
    data = {}
    for i, entry in enumerate(selected_compounds):
        entered_value = entry.get()
        data[f"Item_{i+1}"] = entered_value  # Using a generic key name, you may adjust as needed
        print(entered_value)

# Read items from Excel spreadsheet
excel_file_path = 'C:/Users/4un/Code/TC_Comparison_Code/TC_Compound_Data.xlsx'  # Replace with the path to your Excel file
df = pd.read_excel(excel_file_path)

# Store the data in a matrix (2D list)
matrix = df.values.tolist()

# Create a dictionary to store matrix elements based on header and first column value
matrix_dict = {(header, first_col_value): matrix[row][col]
               for row, first_col_value in enumerate(df.iloc[:, 0])
               for col, header in enumerate(df.columns)}

# Create the main window
root = tk.Tk()
root.title("Thermal Conductivity Model Visualizer")

# Create variables to store checkbox states
compounds = df.iloc[:, 0].dropna().tolist()
checkbox_vars = [IntVar(value=0) for _ in range(len(compounds))]
selected_compounds = []  # List to store selected compounds

# Create checkboxes and add them to the window
for i, item in enumerate(compounds):
    checkbox = tk.Checkbutton(root, text=item, variable=checkbox_vars[i])
    checkbox.pack(anchor=tk.W)
    checkbox.deselect()  # Manually deselect the checkboxes

# Create a button to show selected items
show_button = tk.Button(root, text="Input Composition", command=show_selected)
show_button.pack(pady=10)

# Create a frame to display selected items with Entry widgets
display_frame = tk.Frame(root)
display_frame.pack(pady=10)

# You can adjust the coordinates based on your specific setup
root.geometry("+1920+0")  # For a setup with two 1920x1080 monitors, placing it on the right monitor

# Create a button to save entered data
save_button = tk.Button(root, text="Select Models", command=save_data)
save_button.pack(pady=10)

# Start the main loop
root.mainloop()
