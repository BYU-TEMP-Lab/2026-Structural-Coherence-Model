import numpy as np
from numpy import array
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
import xlrd
import xlwt
from tkinter import *
from tkinter import ttk, filedialog
import pandas as pd
from TC_models import *
import random
import seaborn as sns
import re
import json
from tkinter import filedialog, END
import os
from datetime import datetime


# Global variables to store compound composition and temperature range
composition_inputs = []
temp_inputs = []
selected_items = []
experimental_data = []
measurement_data = []
radio_select = [0]

def save_settings():
    try:
        settings = {
            'selected_compounds': [compounds[i] for i, var in enumerate(checkbox_vars) if var.get()],
            'temp_range': {
                'min_temp': TempMin_entry.get(),
                'max_temp': TempMax_entry.get()
            },
            'selected_methods': [method for method, var in method_vars.items() if var.get() == 1],
            'selected_MSTDB': {key: var.get() for key, var in choices.items() if var.get()},
            'selected_measurement_data': {key: var.get() for key, var in choices_Measurement.items() if var.get()},
            'use_data': run_style.get(),
            'compositions': [[entry.get() for entry in entries] for entries in composition_inputs]  # Ensure proper extraction
        }

        save_path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=(("JSON files", "*.json"), ("All files", "*.*")))
        if save_path:
            with open(save_path, 'w') as f:
                json.dump(settings, f)
            print(f"Settings saved to {save_path}")
    except Exception as e:
        print(f"Error saving settings: {e}")



def load_settings():

    load_path = filedialog.askopenfilename(filetypes=(("JSON files", "*.json"), ("All files", "*.*")))
    if load_path:
        with open(load_path, 'r') as f:
            settings = json.load(f)

        # Load selected compounds
        for i, compound in enumerate(compounds):
            checkbox_vars[i].set(compound in settings['selected_compounds'])
        
        # Load temperature range
        TempMin_entry.delete(0, END)
        TempMin_entry.insert(0, settings['temp_range']['min_temp'])
        TempMax_entry.delete(0, END)
        TempMax_entry.insert(0, settings['temp_range']['max_temp'])
        
        # Load selected methods
        for method in method_vars:
            method_vars[method].set(method in settings['selected_methods'])
        
        # Load MSTDB selections
        for key, var in choices.items():
            var.set(settings['selected_MSTDB'].get(key, 0))
        
        # Load measurement data selections
        for key, var in choices_Measurement.items():
            var.set(settings['selected_measurement_data'].get(key, 0))
        
        # Load use data or model exclusively setting
        run_style.set(settings['use_data'])

        # Load compositions
        InputComposition()  # Call to recreate the entry widgets based on the selected compounds
        for entries, saved_compositions in zip(composition_inputs[-1], settings['compositions']):
            pass
            # for entry, value in zip(entries, saved_compositions):
            #     entry.delete(0, END)
            #     entry.insert(0, value)
        
        print(f"Settings loaded from {load_path}")



def InputComposition(selected_order=None):
    global composition_inputs  # Ensure we are modifying the global variable

    # Get selected compounds from checkboxes
    selected_indices = [i for i, var in enumerate(checkbox_vars) if var.get()]
    
    # Use the provided order if available, otherwise use the checkbox order
    if selected_order is not None:
        selected_items = selected_order
    else:
        selected_items = [compounds[i] for i in selected_indices]

    # Clear the frame before updating with new selections
    for widget in composition_frame.winfo_children():
        widget.destroy()

    # Create Entry widgets only for the selected items
    entry_widgets = []
    for i, item in enumerate(selected_items):
        entry = ttk.Entry(composition_frame, width=10)
        label = ttk.Label(composition_frame, text=f"{item}")
        entry.grid(row=i, column=0, padx=5, pady=2)
        label.grid(row=i, column=1, pady=2)
        entry_widgets.append(entry)

    # Store the entry widgets in the selected_compounds list
    composition_inputs = [entry_widgets]  # Update the global composition_inputs list


# Button function to select a file with experimental data (not used yet)
def SaveData():
    filename = filedialog.askopenfilename(filetypes=(("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls"), ("All files", "*.*")))
    if filename:
        # Do something with the selected file, such as displaying its path or loading its contents
        print("Selected file:", filename)

def handle_radio_selection():
    selected_option = run_style.get()

    if selected_option == 'DATA':
        use_data = 1
    elif selected_option == 'MODEL':
        use_data = 0
    
    radio_select.append(use_data)

def StoreSelectedMSTDB(choices):

    experimental_data.clear()
    for name, var in choices.items():
        if var.get() == 1:
            rows = (TC_only_df[TC_only_df['Formula'] == name])                #print("%s: %s" % (name, var.get()))
            for i, row in rows.iterrows():
                added_row = [(row.iloc[0]),(row.iloc[2]),(row.iloc[1])]
                experimental_data.append(added_row)

    # Keep the menu open after selection
    menu.post(menubutton.winfo_rootx(), menubutton.winfo_rooty() + menubutton.winfo_height())

def StoreSelectedMeasurementData(choices_Measurement):

    measurement_data.clear()
    for name, var in choices_Measurement.items():
        if var.get() == 1:
            measurement_data.append(name)
            # rows = (TC_Measurement_df[TC_Measurement_df['Source'] == name])      
            # for i, row in rows.iterrows():
            #     added_row = [(row.iloc[1])]
            #     measurement_data.append(added_row)

    # Keep the menu open after selection
    menuExp.post(menubuttonExp.winfo_rootx(), menubuttonExp.winfo_rooty() + menubuttonExp.winfo_height())

def ExperimentalData(AB,T):
    
    exp_TC = np.zeros(len(T))
    for j in range(len(T)):
        exp_TC[j] = AB[0] + AB[1]*T[j]

    return(T,exp_TC)


def select_scl_composition(event):
    selected_index = scl_listbox.curselection()
    if selected_index:
        selected_composition = SCL_compositions[selected_index[0]]
        
        # Remove the source from the composition string if present (e.g., '1.0MgCl2 (Lu, 2021)' -> '1.0MgCl2')
        composition_str = selected_composition.split(' (', 1)[0]
        compounds = composition_str.split('-')
        
        # Create lists to store the compound names and values in the order they appear
        compound_order = []
        compound_values = {}
        
        # Parse the compound string (e.g., '0.38MgCl2' -> name='MgCl2', value='0.38')
        for compound in compounds:
            match = re.match(r'(\d+\.?\d*)([A-Za-z0-9]+)', compound.strip())
            if match:
                value = match.group(1)
                name = match.group(2)
                compound_order.append(name)
                compound_values[name] = value
        
        # Get the list of all available compounds from the dataframe
        all_compounds = TC_C_df.iloc[:, 0].dropna().tolist()
        
        # Clear current selections
        for var in checkbox_vars:
            var.set(0)
        
        # Select checkboxes for all compounds in the composition
        for name in compound_order:
            if name in all_compounds:
                index = all_compounds.index(name)
                checkbox_vars[index].set(1)
        
        # Trigger composition input update with the selected order
        InputComposition(selected_order=compound_order)
        
        # Set composition values in the correct order as percentages
        if composition_inputs and composition_inputs[-1]:
            for i, entry in enumerate(composition_inputs[-1]):
                if i < len(compound_order):
                    compound_name = compound_order[i]
                    value = float(compound_values.get(compound_name, 0)) * 100  # Convert to percentage
                    entry.delete(0, END)
                    entry.insert(0, f"{value:.2f}".rstrip('0').rstrip('.'))  # Format to remove trailing zeros



# Button function to run the estimation model codes and graph everything (plan to develop more options)
def CreateGraphs():
    # Get selected compounds from checkboxes
    selected_indices = [i for i, var in enumerate(checkbox_vars) if var.get()]
    selected_items = [compounds[i] for i in selected_indices]

    # Store the temp range entry widgets in a list
    temp_inputs.extend([TempMin_entry, TempMax_entry])

    # Store the selected methods in a list
    selected_methods = [method for method, var in method_vars.items() if var.get() == 1]

    # Retrieve and process composition inputs, ensuring they match the selected compounds order
    # Create a list of (compound, value) pairs to maintain the correct association
    compound_pairs = list(zip(selected_items, [float(entry.get())*0.01 for entry in composition_inputs[-1]]))
    
    # Sort the pairs to ensure consistent ordering (alphabetical by compound name)
    compound_pairs.sort(key=lambda x: x[0])
    
    # Extract the sorted compounds and values
    selected_items_sorted = [pair[0] for pair in compound_pairs]
    compound_values = [pair[1] for pair in compound_pairs]

    # Retrieve temperature range data from widgets
    min_temp = float(temp_inputs[0].get())
    max_temp = float(temp_inputs[1].get())
    temp_range = [min_temp,max_temp]

    # Apply publication-style matplotlib settings
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    # Legend and small text/annotation sizing (~0.85x of tick labels)
    _tick_size = 12
    _small_text = int(round(0.85 * _tick_size))  # 10 when tick labels are 12
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.fontsize'] = _small_text
    
    GraphFig = plt.figure()
    ax = GraphFig.add_subplot(111)
    # Thicker spines to match publication style
    for side in ['left', 'bottom', 'right', 'top']:
        try:
            ax.spines[side].set_linewidth(1.5)
        except Exception:
            pass
    plt.xlabel('Temperature (K)')
    # Use middle dots and superscripts via mathtext
    plt.ylabel('Thermal Conductivity (W·m$^{-1}$·K$^{-1}$)')

    ModelMixtureLabel = []
    ModelLabels = []
    DataLabels = []

    # Define a list of shape markers
    shape_markers = ['o', 's', '^', 'x', '+', 'D', 'v', '*', 'P', '>', '<', 'h', 'H', 'd', '|', '_', '8', '.', ',', '1', '2', '3', '4', 'p', 'X', 'o', 's', '^', 'x', '+', 'D', 'v', '*', 'P', '>', '<', 'h', 'H', 'd', '|', '_', '8', '.', ',', '1', '2', '3', '4', 'p', 'X']
    
    # Define the number of colors needed (assuming len(measurement_data) is the number of iterations)
    num_colors = len(selected_methods)+len(experimental_data)+len(measurement_data)

    palette = sns.color_palette('deep', num_colors)

    color_i = 0
    # Dictionary to store model predictions at melting temperature
    model_predictions_at_melt = {}
    # Dictionary to store model predictions at minimum measured temperature
    model_predictions_at_min_temp = {}
    # Dictionary to store model results to avoid recalculation
    model_results = {}
    # Variable to store minimum measured temperature
    min_measured_temp = None
    
    # Dictionary to store additional model outputs (specific heat, sound velocity, etc.)
    additional_outputs = {}
    
    for i, method in enumerate(selected_methods):
        run = estimation_methods_dict[method]
        # Run the model and store results
        if radio_select[-1] == 0:
            # Use the sorted compounds and values for model calculation
            result = run(TC_C_df, MSTDB_df, SCL_PDF_df, selected_items_sorted, compound_values, temp_range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0)
        else:
            # Use the sorted compounds and values for model calculation
            result = run(TC_C_df, MSTDB_df, SCL_PDF_df, selected_items_sorted, compound_values, temp_range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=1)
        
        # Check if the result is a tuple (T, dict) or the old format (T, array)
        if isinstance(result, tuple) and len(result) == 2 and isinstance(result[1], dict):
            T, result_dict = result
            lambda_mix_T = result_dict['thermal_conductivity']
            
            # Store additional outputs if they exist
            if method in ['Present Model', 'Present Model, Mix Data']:
                additional_outputs[f'{method}_specific_heat_m'] = result_dict.get('specific_heat_m', '')
                additional_outputs[f'{method}_specific_heat_prime'] = result_dict.get('specific_heat_prime', '')
                additional_outputs[f'{method}_sound_velocity_m'] = result_dict.get('sound_velocity_m', '')
                additional_outputs[f'{method}_sound_velocity_prime'] = result_dict.get('sound_velocity_prime', '')
        else:
            # Handle old format for backward compatibility
            T, lambda_mix_T = result
        
        # Store the full model results for this method
        model_results[method] = (T.copy(), lambda_mix_T.copy())

        # Store the prediction at melting temperature for CSV saving
        melting_temp = float(temp_inputs[0].get()) if len(temp_inputs) > 0 else min_temp
        # Find the closest temperature index to melting temperature
        if len(T) > 0 and len(lambda_mix_T) > 0:
            closest_idx = np.argmin(np.abs(T - melting_temp))
            if isinstance(lambda_mix_T, dict):
                model_predictions_at_melt[method] = float(lambda_mix_T['thermal_conductivity'][closest_idx])
            else:
                model_predictions_at_melt[method] = float(lambda_mix_T[closest_idx])
        else:
            model_predictions_at_melt[method] = ''
            
        # Store the prediction at minimum measured temperature (will be set later)
        model_predictions_at_min_temp[method] = ''

        # Format formula label using the sorted compounds and values
        label = '-'.join([f"{comp}{name}" for comp, name in zip(compound_values, selected_items_sorted)])
        formatted_label = f"{label}, {method}"

        # Generate a new color for each function
        color = palette[color_i]
        color_i = color_i + 1
        
        # Plot the thermal conductivity
        if isinstance(lambda_mix_T, dict):
            plt.plot(T, lambda_mix_T['thermal_conductivity'], label=formatted_label, color=color)
        else:
            plt.plot(T, lambda_mix_T, label=formatted_label, color=color)
            
        plt.tight_layout()

        if label not in ModelMixtureLabel:
            ModelMixtureLabel = np.append(ModelMixtureLabel,label)

        if method not in ModelLabels:
            ModelLabels = np.append(ModelLabels,method)

    if len(experimental_data) == 0:
        #print("No MSTDB data selected")
        pass
    else:
        T = np.linspace(temp_range[0],temp_range[1],num=100)
        for i in range(len(experimental_data)):
            T,exp_TC = ExperimentalData(experimental_data[i][1],T)
            # Generate a new color for each function
            color = palette[color_i]
            color_i = color_i+1       
            plt.plot(T,exp_TC,label=((experimental_data[i][0]),(experimental_data[i][2])),color=color,linestyle='--')
            plt.tight_layout()
            
            if experimental_data[i][0] not in DataLabels:
                DataLabels = np.append(DataLabels,experimental_data[i][0])
        print(experimental_data[i][0], " MSTDB avg. at melt temp: ",exp_TC[0])
    
    if len(measurement_data) == 0:
        pass#print("No measurement data selected")
    else:
        AvgMeltLambdaList = []
        # Store experimental values at melting temperature for CSV saving
        experimental_at_melt = []
        # List to store all minimum temperatures from measurement data
        all_min_temps = []
        for i in range(len(measurement_data)):
            rows = TC_Measurement_df[TC_Measurement_df['Source'] == measurement_data[i]]
            T = np.zeros(len(rows))
            meas_TC = np.zeros(len(rows))
            for j, (_, row) in enumerate(rows.iterrows()):
                T[j],meas_TC[j] = [(row.iloc[3]),(row.iloc[4])]
            # Store the minimum temperature from this measurement set
            if len(T) > 0:
                all_min_temps.append(min(T))
            # Prepare linear fit
            coefficients = np.polyfit(T,meas_TC,1)
            linear_fit = np.poly1d(coefficients)
            x_fit = np.linspace(min(T),max(T),100)
            #x_fit = np.linspace(min_temp,max(T),100)
            # Get variance
            std = np.std(meas_TC)
            # Generate a new color for each function
            color = palette[color_i]
            color_i = color_i+1    
            ###  
            # if 'Merritt' in measurement_data[i]:
            #     shape_markers[i] = 'o'
            #     plt.scatter(T,meas_TC,label=(measurement_data[i]),color='blue',marker=shape_markers[i])
            # else:
            #     shape_markers[i] = '+'
            #     plt.scatter(T,meas_TC,label=(measurement_data[i]),color=color,marker=shape_markers[i])
            # if 'KCl' in measurement_data[i]:
            #     shape_markers[i] = '^'
            #     plt.scatter(T,meas_TC,label=(measurement_data[i]),color=color,marker=shape_markers[i])
            # elif 'NaF':
            #     shape_markers[i] = 'x'
            #     plt.scatter(T,meas_TC,label=(measurement_data[i]),color=color,marker=shape_markers[i])
            # else:
            #     shape_markers[i] = '^'
            #     plt.scatter(T,meas_TC,label=(measurement_data[i]),color=color,marker=shape_markers[i])
            # elif 'NaCl':
            #     shape_markers[i] = 'x'
            #     pattern = r'\((.*?)\)'
            #     match = re.search(pattern, measurement_data[i])
            #     if match:
            #         label = match.group(1)
            #     plt.scatter(T,meas_TC,label=f"FLiNaK ({label})",color=color,marker=shape_markers[i])
            # else:
            #     shape_markers[i] = '^'
            ###
            
            # Calculate TC at melting temperature instead of minimum temperature
            melting_temp = float(temp_inputs[0].get()) if len(temp_inputs) > 0 else min_temp
            tc_at_melt = linear_fit(melting_temp)
            AvgMeltLambdaList.append(tc_at_melt)
            experimental_at_melt.append(tc_at_melt)
            plt.scatter(T,meas_TC,label=(measurement_data[i]),color=color,marker=shape_markers[i])
            plt.plot(x_fit,linear_fit(x_fit),color=color,linestyle='dotted')
            plt.fill_between(x_fit, linear_fit(x_fit)-std, linear_fit(x_fit)+std, alpha=0.12,color=color)
            #print("StdDev, ", measurement_data[i], " : ",std)
            
            plt.tight_layout() 

            Data = measurement_data[i].split()[0]
            if Data not in DataLabels:
                DataLabels = np.append(DataLabels,Data)

        AvgMeltLambda = np.mean(AvgMeltLambdaList)
        print(measurement_data[i][0]," Exp avg. at melt temp: ", AvgMeltLambda)
        
        # Set the global minimum measured temperature (minimum of all measurement sets)
        if len(all_min_temps) > 0:
            min_measured_temp = min(all_min_temps)
            
            # Calculate model predictions at minimum measured temperature using stored results
            for method in selected_methods:
                if method in model_results:
                    T_min, lambda_mix_T_min = model_results[method]
                else:
                    # Fallback to running the model if results aren't stored (shouldn't happen)
                    run = estimation_methods_dict[method]
                    if radio_select[-1] == 0:
                        T_min,lambda_mix_T_min = run(TC_C_df,MSTDB_df,SCL_PDF_df,selected_items,compound_values,temp_range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0)
                    else:
                        T_min,lambda_mix_T_min = run(TC_C_df,MSTDB_df,SCL_PDF_df,selected_items,compound_values,temp_range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=1)
                
                # Find the closest temperature index to minimum measured temperature
                if len(T_min) > 0 and len(lambda_mix_T_min) > 0:
                    closest_idx_min = np.argmin(np.abs(T_min - min_measured_temp))
                    model_predictions_at_min_temp[method] = float(lambda_mix_T_min[closest_idx_min])
                else:
                    model_predictions_at_min_temp[method] = ''   

    #plt.ylim(BOTTOM = 0)
    plt.legend()

    # Get the current date
    current_date = datetime.now().strftime("%Y%m%d")
    save_title = f'{label}_{selected_methods}_{current_date}'

    ModelMixtureLabel = [str(x) for x in ModelMixtureLabel]
    ModelLabels = [str(x) for x in ModelLabels]
    DataLabels = [str(x) for x in DataLabels]

    # Create a formatted title
    save_title = f"[{'_'.join(ModelMixtureLabel)}]_[{'_'.join(ModelLabels)}]_[{'_'.join(DataLabels)}]"


    Save_fig(save_title)
    plt.show()
    # Use compound_values and temperature range data as needed
    # print("Compound Data:", (compound_values[0]+compound_values[1]))
    # print("Temperature Range Data:", (temp_range[0]+temp_range[1]))
    # print("Selected Estimation Methods:", selected_methods)
    # --- BEGIN: Collect melt point results for saving ---
    melt_results = []
    # Prepare shared fields for this composition
    composition_str = ','.join(ModelMixtureLabel) if 'ModelMixtureLabel' in locals() else ''
    melting_temp = float(temp_inputs[0].get()) if len(temp_inputs) > 0 else ''
    
    # Get PDF source from SCL data for this composition
    pdf_source = ''
    if composition_str:
        # Get the currently selected SCL composition if any
        selected_scl_index = scl_listbox.curselection()
        if selected_scl_index:
            selected_scl = scl_listbox.get(selected_scl_index[0])
            # Extract the source from the selected SCL composition (e.g., '1.0MgCl2 (Lu, 2021)' -> 'Lu, 2021')
            if '(' in selected_scl and ')' in selected_scl:
                source = selected_scl.split('(', 1)[1].rsplit(')', 1)[0].strip()
                # Look for matching composition and source in SCL_PDF_df
                matching_rows = SCL_PDF_df[(SCL_PDF_df.iloc[:, 0] == composition_str) & 
                                         (SCL_PDF_df.iloc[:, 1] == source)]
                if not matching_rows.empty:
                    pdf_source = source
    # Prepare experimental (measurement) results for this composition
    tc_exp = ''
    exp_ref = ''
    if len(measurement_data) > 0 and 'experimental_at_melt' in locals() and len(experimental_at_melt) > 0:
        exp_ref = ','.join([str(ref) for ref in measurement_data])
        # Compute the average if multiple values (now at melting temperature)
        try:
            vals = [float(val) for val in experimental_at_melt]
            tc_exp = str(sum(vals) / len(vals)) if len(vals) > 0 else ''
        except Exception:
            tc_exp = ''
    # Prepare MSTDB results for this composition
    tc_mstdb = ''
    mstdb_ref = ''
    if len(experimental_data) > 0 and 'exp_TC' in locals() and len(exp_TC) > 0:
        mstdb_ref = ','.join([str(experimental_data[i][2]) if len(experimental_data[i]) > 2 else 'MSTDB' for i in range(len(experimental_data))])
        # Compute the average if multiple values
        try:
            vals = [float(exp_TC[0])]  # extend here if exp_TC can be a list of multiple MSTDB values
            tc_mstdb = str(sum(vals) / len(vals)) if len(vals) > 0 else ''
        except Exception:
            tc_mstdb = ''
    # Build one row for the composition, with a column for each model from functionlibrary()
    import importlib.util
    import sys
    # Dynamically import TC_models.py to get functionlibrary
    model_module_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TC_models.py')
    spec = importlib.util.spec_from_file_location("TC_models", model_module_path)
    model_mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = model_mod
    spec.loader.exec_module(model_mod)
    all_model_names = list(model_mod.functionlibrary().keys())
    # Map model name to its TC value using the stored predictions
    model_results = {name: '' for name in all_model_names}
    # Use the stored model predictions from the execution loop
    if 'model_predictions_at_melt' in locals():
        for method in model_predictions_at_melt:
            if method in model_results:
                model_results[method] = model_predictions_at_melt[method]
    # Construct the row as a dict
    melt_row = {
        'composition': composition_str,
        'pdf_source': pdf_source,
        'melting_temp': melting_temp,
        'min_measured_temp': min_measured_temp if 'min_measured_temp' in locals() and min_measured_temp is not None else '',
        'tc_experimental': tc_exp,
        'exp_reference': exp_ref,
        'tc_mstdb': tc_mstdb,
        'mstdb_reference': mstdb_ref
    }
    # Add model columns at melting temperature
    for model_name in all_model_names:
        melt_row[f"{model_name}_at_melt"] = model_results[model_name]
    # Add model columns at minimum measured temperature
    for model_name in all_model_names:
        if 'model_predictions_at_min_temp' in locals() and model_name in model_predictions_at_min_temp:
            melt_row[f"{model_name}_at_min_temp"] = model_predictions_at_min_temp[model_name]
        else:
            melt_row[f"{model_name}_at_min_temp"] = ''
    
    # Add additional outputs (specific heat, sound velocity, etc.) to the melt_row
    if 'additional_outputs' in locals() and additional_outputs:
        for key, value in additional_outputs.items():
            melt_row[key] = value
    
    melt_results.append(melt_row)
    window.melt_results = melt_results
    # --- END: Collect melt point results for saving ---

    # If the save results checkbox is checked, save to CSV
    if save_results_var.get():
        save_results_to_csv()

def Save_fig(title):
    # Split the title into components
    parts = title.split(']_[')
    
    # Process composition part (first part)
    if parts:
        # Extract and round composition values
        comp_parts = parts[0].replace('[', '').split('_')
        rounded_comps = []
        for comp in comp_parts:
            try:
                # Try to round the number if it's a float
                num = float(comp)
                rounded = round(num, 5)
                # Remove trailing .0 if it's an integer after rounding
                if rounded == int(rounded):
                    rounded = int(rounded)
                rounded_comps.append(str(rounded))
            except ValueError:
                # If not a number, keep as is
                rounded_comps.append(comp)
        parts[0] = '_'.join(rounded_comps)
    
    # Simplify model and data labels (second and third parts)
    if len(parts) > 1:
        # Just take the first label from each group
        parts[1] = parts[1].split('_')[0]  # First model
        if len(parts) > 2:
            parts[2] = parts[2].split('_')[0]  # First data source
    
    # Create a simpler filename
    current_date = datetime.now().strftime("%Y%m%d")
    save_title = f'TC_{parts[0]}_{current_date}.png'
    
    # Ensure the filename isn't too long (max 200 chars)
    if len(save_title) > 200:
        save_title = save_title[:200] + '.png'
    
    folder_name = "Comparison_Plot_Figs"
    os.makedirs(os.path.join(os.getcwd(), folder_name), exist_ok=True)
    save_path = os.path.join(os.getcwd(), folder_name, save_title)

    # Save the plot
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Figure saved to: {save_path}")
    return

# TC_composition_inputs - Convert excel and CSV files
# Read the compound data, handling potential source information in the compound names
TC_C_df = pd.read_excel('TC_compound_data.xlsx')
# Ensure the first column is treated as string to preserve any source information
TC_C_df.iloc[:, 0] = TC_C_df.iloc[:, 0].astype(str)

# MSTDB-TP_inputs - Convert excel and CSV files
MSTDB_df = pd.read_csv('MSTDB.csv')

# TC_measurement_data inputs
TC_Measurement_df = pd.read_excel('TC_measurement_data.xlsx')

# SCL_PDF_Analyze results to input
SCL_PDF_df = pd.read_csv('scl_results.csv')

# Create background window
window = Tk()
window.title("Thermal Conductivity Estimation Method Comparison")
window.geometry('580x600')
window.columnconfigure(0,weight=1)
window.rowconfigure(0,weight=1)

# Create main GUI window
window_frame = ttk.Frame(window)
window_frame.pack(fill='both', expand=True) 

# Create left, middle, and right portions
left_frame = ttk.Frame(window_frame, borderwidth=2, relief="solid")
left_frame.pack(side='left', expand=True, fill='both')
middle_frame = ttk.Frame(window_frame, borderwidth=2, relief="solid")
middle_frame.pack(side='left', expand=True, fill='both')
right_frame = ttk.Frame(window_frame, borderwidth=2, relief="solid")
right_frame.pack(side='left', expand=True, fill='both')

# Create label for compound selection
ttk.Label(left_frame,text="Select Compounds").pack(anchor=N,pady=3)


# Create variables to store checkbox states for compounds
# Get all compounds, ensuring we handle any source information in the names
compounds = TC_C_df.iloc[:, 0].dropna().tolist()
# Create a mapping of display names to original names to preserve source information
compound_display_names = []
for compound in compounds:
    # If the compound has a source in parentheses, use it as is
    if '(' in compound and ')' in compound:
        display_name = compound
    # Otherwise, check if there's a source column and append it
    elif len(TC_C_df.columns) > 1 and 'Source' in TC_C_df.columns:
        source = TC_C_df.loc[TC_C_df.iloc[:, 0] == compound, 'Source'].iloc[0] if not pd.isna(TC_C_df.loc[TC_C_df.iloc[:, 0] == compound, 'Source']).any() else ''
        display_name = f"{compound} ({source})" if source else compound
    else:
        display_name = compound
    compound_display_names.append(display_name)

checkbox_vars = [IntVar(value=0) for _ in range(len(compounds))]

# Create checkboxes and add them to the window
for i, item in enumerate(compounds):
    checkbox = ttk.Checkbutton(left_frame, text=item, variable=checkbox_vars[i]).pack(anchor=W,padx=5)

# Create vertical frame for temperature range selection
temp_frame = ttk.Frame(middle_frame)
temp_frame.pack(expand=True, fill=X)
ttk.Label(temp_frame,text="Select Temp Range (K)").pack(anchor=S,pady=3)

# Create horizontal frame for temp min and max entries and labels
temp_entry_frame = ttk.Frame(temp_frame)
temp_entry_frame.pack(expand=True)

# Create vertical frame for temp min entry box and label
min_temp_entry_frame = ttk.Frame(temp_entry_frame)
min_temp_entry_frame.pack(side='left',expand=True, padx=2, pady=2)

# Create empty textbox for min temp  input
TempMin = float
TempMin_entry = ttk.Entry(min_temp_entry_frame,width=10,textvariable=TempMin)
TempMin_entry.pack()
ttk.Label(min_temp_entry_frame,text="Min Temp").pack()

# Create vertical frame for temp max entry box and label
max_temp_entry_frame = ttk.Frame(temp_entry_frame)
max_temp_entry_frame.pack(side='left',expand=True, padx=2, pady=2)

# Create empty textbox for max temp  input
TempMax = float
TempMax_entry = ttk.Entry(max_temp_entry_frame,width=10,textvariable=TempMax)
TempMax_entry.pack()
ttk.Label(max_temp_entry_frame,text="Max Temp").pack()

# Create variables to store estimation method checkbox states
estimation_methods_dict = functionlibrary()
method_vars = {method: IntVar() for method in estimation_methods_dict}

#Create checkboxes for estimation methods
method_frame = ttk.Frame(middle_frame)
method_frame.pack(expand=True, fill='both')
ttk.Label(method_frame,text="Select Methods").pack(pady=3)
for i, (method, description) in enumerate(estimation_methods_dict.items()):
    checkbox = ttk.Checkbutton(method_frame, text=f"{method}", variable=method_vars[method]).pack(anchor=W,padx=5)

# Create a button to input composition for selected compounds
ttk.Button(middle_frame, text="Input Composition (mol%)", command=InputComposition).pack(anchor=S,pady=3)

# Create a frame to display selected compounds with empty textboxes
composition_frame = ttk.Frame(middle_frame)
composition_frame.pack(expand=True, fill='both')

### Right Frame ##################################################################

# Extract all SCL compositions including their sources
SCL_compositions = []
for idx, row in SCL_PDF_df.iterrows():
    comp_name = row.iloc[0]
    source = row.iloc[1] if len(row) > 1 else 'Unknown'
    # Append source to composition name if it's not already included
    if not comp_name.endswith(f' ({source})'):
        comp_name = f'{comp_name} ({source})'
    SCL_compositions.append((comp_name, idx))  # Store both display name and original index

# Sort the compositions by their original row index to maintain order
SCL_compositions = [comp[0] for comp in sorted(SCL_compositions, key=lambda x: x[1])]

# Create a frame for SCL compositions
scl_frame = ttk.Frame(right_frame, borderwidth=2, relief="solid")
scl_frame.pack(side='top', expand=True, fill='both', pady=10)

ttk.Label(scl_frame, text="SCL Compositions").pack(pady=5)

# Create a listbox to display SCL compositions
scl_listbox = tk.Listbox(scl_frame, width=20, height=5)
scl_listbox.pack(expand=True, fill='both', padx=5, pady=5)

# Populate the listbox with unique SCL compositions (only most recent)
for composition in SCL_compositions:
    scl_listbox.insert(tk.END, composition)

# Add a scrollbar to the listbox
scl_scrollbar = ttk.Scrollbar(scl_frame, orient="vertical", command=scl_listbox.yview)
scl_scrollbar.pack(side="right", fill="y")
scl_listbox.configure(yscrollcommand=scl_scrollbar.set)

# Bind the selection event to the listbox
scl_listbox.bind('<<ListboxSelect>>', select_scl_composition)

# Create drop-down list of MSTDB-TP data to select for comparison
menubutton = tk.Menubutton(right_frame, text="MSTDB-TP Data", 
                           indicatoron=True, borderwidth=1, relief="raised")
menu = tk.Menu(menubutton, tearoff=False)
menubutton.configure(menu=menu)
menubutton.pack(padx=10, pady=15)

# Initialize lists to store extracted data
salt = []
compositions = []
thermal_conductivities = []
temps = []
references = []
melting_T = []

# Iterate over each row in the DataFrame
for index, row in MSTDB_df.iterrows():
    # Check if Thermal Conductivity appears to be a number
    conductivity_A = str(row['Thermal Conductivity (W/m K):  A + B*T(K)'])
    conductivity_B = str(row.iloc[row.index.get_loc('Thermal Conductivity (W/m K):  A + B*T(K)') + 1])
    try:
        conductivityA = pd.to_numeric(conductivity_A, errors='coerce')
        conductivityB = pd.to_numeric(conductivity_B, errors='coerce')
        if not pd.isnull(conductivityA):
            # Split compound names and compositions
            df_compunds = str(row['Formula']).split('-')
            temp_range = str(row.iloc[row.index.get_loc('Thermal Conductivity (W/m K):  A + B*T(K)') + 2]).split('-')
            reference = str(row.iloc[row.index.get_loc('Thermal Conductivity (W/m K):  A + B*T(K)') + 4])

            composition_str = str(row['Composition (Mole %)'])

            meltingT = str(row['Melting T (K)'])
            
            # Check if it's a single compound
            if '-' not in composition_str:
                composition_percentages = [1.0]  # Pure salt, 100%
                salt.append(df_compunds[0])  # Single compound, no need for coefficients
            else:
                composition_percentages = list(map(float, composition_str.split('-')))
                salt.append('-'.join([f"{comp}{salt}" for comp, salt in zip(composition_percentages, df_compunds)]))
          
            if not pd.isnull(conductivityB):
                thermal_conductivities.append([float(conductivity_A), float(conductivity_B)]) 
            else:
                thermal_conductivities.append([float(conductivity_A), 0])

            compositions.append(composition_percentages)
            temps.append(temp_range) 
            references.append(reference)
            melting_T.append(meltingT)         
                    
    except ValueError:
        pass
        
# Create a new DataFrame with the extracted data
TC_only_df = pd.DataFrame({'Formula': salt, 'Reference': references, 'Thermal Conductivity': thermal_conductivities, 'Melting T (K)': melting_T})

choices = {}

# Loop through each row in the DataFrame
for index, row in TC_only_df.iterrows():
    # Extract the 'Formula' and 'Melting T (K)' values from the row
    formula = row['Formula']
    melting_t = row['Melting T (K)']
    
    # Create a Tkinter IntVar and add it to the choices dictionary
    choices[formula] = tk.IntVar(value=0)
    
    # Add a checkbutton to the menu with the formula and melting_t as the label
    menu.add_checkbutton(label=f"{formula}, ({melting_t} K)",
                         variable=choices[formula],
                         onvalue=1, offvalue=0,
                         command=lambda: StoreSelectedMSTDB(choices))

    
    
# Create drop-down list of experimental data to select for comparison
menubuttonExp = tk.Menubutton(right_frame, text="Experimental Data", 
                           indicatoron=True, borderwidth=1, relief="raised")
menuExp = tk.Menu(menubuttonExp, tearoff=False)
menubuttonExp.configure(menu=menuExp)
menubuttonExp.pack(padx=10, pady=15)

choices_Measurement = {}
options = TC_Measurement_df[['Dataset ID', 'Source',]].drop_duplicates(subset='Dataset ID').iloc[:,1].values
for choice in options:
    choices_Measurement[choice] = tk.IntVar(value=0)
    menuExp.add_checkbutton(label=choice, variable=choices_Measurement[choice], 
                         onvalue=1, offvalue=0, 
                         command=lambda: StoreSelectedMeasurementData(choices_Measurement))
    

# Create option to use available data or run exlusively with models
radio_frame = ttk.Frame(right_frame)
radio_frame.pack(anchor=N,expand=True,fill=X,pady=4)
run_style = StringVar()
home = ttk.Radiobutton(radio_frame, text='Use available data', variable=run_style, value='DATA', command=handle_radio_selection).pack(anchor=W,padx=5)
office = ttk.Radiobutton(radio_frame, text='Use models exlusively', variable=run_style, value='MODEL', command=handle_radio_selection).pack(anchor=W,padx=5)


# # Create option to save data
# fileexplorer_icon = tk.PhotoImage(file="fileexplorer_icon.png")
# ttk.Button(right_frame, text="Select Save Location", compound=tk.LEFT, command=SaveData, image=fileexplorer_icon).pack(anchor=S, pady=10)

# Add these buttons to the right_frame or where appropriate in your GUI layout

ttk.Button(right_frame, text="Save Settings", command=save_settings).pack(anchor=S, pady=5)
ttk.Button(right_frame, text="Load Settings", command=load_settings).pack(anchor=S, pady=5)

def save_results_to_csv():
    import csv
    import os
    # Define the path to the CSV in the SAME directory as this script
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TC_calc_results.csv")
    # Always use all model columns from functionlibrary in TC_models.py
    import importlib.util
    import sys
    model_module_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TC_models.py')
    spec = importlib.util.spec_from_file_location("TC_models", model_module_path)
    model_mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = model_mod
    spec.loader.exec_module(model_mod)
    model_columns = list(model_mod.functionlibrary().keys())
    # Define header
    # Create model column names for melting temperature
    model_columns_at_melt = [f"{col} at Melt Temp (W/m-K)" for col in model_columns]
    # Create model column names for minimum measured temperature
    model_columns_at_min = [f"{col} at Min Meas Temp (W/m-K)" for col in model_columns]
    
    # Add specific heat and sound velocity columns for GECM and GECM_Mix
    specific_heat_columns = [
        'GECM c_m (J/kg-K)',
        'GECM c\' (J/kg-K²)',
        'GECM_Mix c_m (J/kg-K)',
        'GECM_Mix c\' (J/kg-K²)'
    ]
    
    sound_velocity_columns = [
        'GECM v_m (m/s)',
        'GECM v\' (m/s/K)',
        'GECM_Mix v_m (m/s)',
        'GECM_Mix v\' (m/s/K)'
    ]
    
    header = [
        'Salt Composition',
        'PDF Source',
        'Melting Temp (K)',
        'Min Measured Temp (K)',
    ] + model_columns_at_melt + model_columns_at_min + [
        'TC Experimental (W/m-K)',
        'Exp. Reference',
        'TC MSTDB (W/m-K)',
        'MSTDB Reference'
    ] + specific_heat_columns + sound_velocity_columns
    # Gather all results to save (models and experiment)
    rows = []
    for result in getattr(window, 'melt_results', []):
        # Get GECM and GECM_Mix specific heat and sound velocity data
        # Extract values directly from the result dictionary
        gecm_c_m = result.get('Present Model_specific_heat_m', '')
        gecm_c_prime = result.get('Present Model_specific_heat_prime', '')
        gecm_mix_c_m = result.get('Present Model, Mix Data_specific_heat_m', '')
        gecm_mix_c_prime = result.get('Present Model, Mix Data_specific_heat_prime', '')
        
        gecm_v_m = result.get('Present Model_sound_velocity_m', '')
        gecm_v_prime = result.get('Present Model_sound_velocity_prime', '')
        gecm_mix_v_m = result.get('Present Model, Mix Data_sound_velocity_m', '')
        gecm_mix_v_prime = result.get('Present Model, Mix Data_sound_velocity_prime', '')
        
        row = [
            result.get('composition', ''),
            result.get('pdf_source', ''),
            result.get('melting_temp', ''),
            result.get('min_measured_temp', ''),
        ]
        # Add model columns at melting temperature
        for model in model_columns:
            row.append(result.get(f"{model}_at_melt", ''))
        # Add model columns at minimum measured temperature
        for model in model_columns:
            row.append(result.get(f"{model}_at_min_temp", ''))
        row += [
            result.get('tc_experimental', ''),
            result.get('exp_reference', ''),
            result.get('tc_mstdb', ''),
            result.get('mstdb_reference', '')
        ]
        
        # Add specific heat data
        row.extend([gecm_c_m, gecm_c_prime, gecm_mix_c_m, gecm_mix_c_prime])
        # Add sound velocity data
        row.extend([gecm_v_m, gecm_v_prime, gecm_mix_v_m, gecm_mix_v_prime])
        rows.append(row)
    if not rows:
        print("No melt point results available to save.")
        return
    # Write or append to CSV
    write_header = not os.path.exists(csv_path)
    with open(csv_path, 'a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(header)
        writer.writerows(rows)
    print(f"Results saved to {csv_path}")

# Add the Save Results checkbox
save_results_var = tk.BooleanVar(value=False)
save_results_chk = ttk.Checkbutton(right_frame, text="Save Results to CSV", variable=save_results_var)
save_results_chk.pack(anchor=S, pady=5)

# Run button
ttk.Button(right_frame, text="Create Graphs", command=CreateGraphs).pack(anchor=S, pady=30)

window.mainloop()


