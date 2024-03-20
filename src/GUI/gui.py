##################
# TODO:
# 1. Separate GUI into processing and plotting sections
# 2. Add plot_conditions
# 3. Add normalization_method
# 4. Add xticks, yticks
# 5. Add color label
# 6. Add plot size
# 7. Add dose concs
# 8. change to loops
##################

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from table import *
import json
import os

HOME_DIR = os.getcwd()
IMAGES_DIR = ""
json_dict = {}
conditions_list = []
condition_names = [""]

# toggle_checkbox shows folder selection widgets if checkbox is checked
def toggle_checkbox():
    if var.get():
        folder_select_lbl.pack(side=TOP, anchor = "w", padx = 10, before = sample_name_label)
        folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10, before = sample_name_label)
        dir_entry.pack(side=TOP, anchor = "w", padx = 10, before = sample_name_label)
    else:
        folder_select_lbl.pack_forget()
        folder_select_lbl_btn.pack_forget()
        dir_entry.pack_forget()

# saves directory to variable, modifies entry widget to show directory string
def images_directory():
    IMAGES_DIR = filedialog.askdirectory()
    dir_entry.insert(0, string=IMAGES_DIR)

# read_sample_name grabs the entry from the sample_name_entry widget
# if the enter_cells widget is inactive (widget to save cells in table)
# then it re-enables the widget, also re-enables the plate_count_option_menu widget
def read_sample_name():
    condition_names.clear()
    sample_name = sample_name_entry.get()
    sample_prompt_text.config(text = f"Enter cells for sample: {sample_name}")
    enter_cells.configure(state="active")
    plate_count_option_menu.configure(state="active")

# save_sample_cells() saves the current list of sample_cells to
# the conditions_list dictionary, then clears the selected cells and the table
def save_sample_cells():
    conditions = {}
    sample_name = sample_name_entry.get()
    conditions[sample_name] = sorted(list(sample_cells))
    conditions_list.append(conditions)
    print(conditions_list)
    if plate_cells_var.get() < plate_count_var.get():
        plate_cells_var.set(plate_cells_var.get() + 1)
    else:
        enter_cells.configure(state="disabled")
    sample_cells.clear()
    table.clear_color()

# disable_plate_count disables the plate_count_option_menu widget, prevents user from editing it after hitting "Enter"
def disable_plate_count():
    plate_count_option_menu.configure(state="disabled")
    plate_cells_var.set(1)

# save all parameters to json file
def save_to_json():
    json_dict["notes"] = notes_entry.get("1.0", "end-1c")
    json_dict["aquisition_frequency"] = int(acquisition_freq.get("1.0", "end-1c"))
    json_dict["directory"] = HOME_DIR
    json_dict["images_directory"] = IMAGES_DIR
    json_dict["conditions"] = conditions_list
    json_dict["sig"] = 2
    json_dict["blockDiameter"] = [501, 101]
    json_dict["shift_thresh"] = 50
    with open(HOME_DIR + "/test_configs.json", "w") as file:
        json.dump(json_dict, file, indent = 4)

def plate_count_select(value):
    plate_count_var.set(value)

if __name__ == "__main__":
    root = Tk()
    root.title("Scanner")
    root.geometry("1000x950")
    frm = ttk.Frame(root)
    frm.grid()
    frm.pack(padx=2, pady=0)

    # EXPERIMENT NOTES ENTRY
    notes_lbl = ttk.Label(root, text = "Experiment Notes")
    notes_lbl.pack(side=TOP, anchor = "w", padx = 10)

    notes_entry = Text(root, width = 35, height = 4)
    notes_entry.pack(side=TOP, anchor = "w", padx = 10)

    # ACQUISITION FREQUENCY
    acquisition_lbl = ttk.Label(root, text = "Acquisition Frequency (#/hr)")
    acquisition_lbl.pack(side=TOP, anchor = "w", padx = 10)

    acquisition_freq = Text(root, width = 35, height = 4)
    acquisition_freq.pack(side=TOP, anchor = "w", padx = 10)

    # CHECKBOX FOR IMAGES DIRECTORY
    var = ttk.IntVar()
    checkbox = ttk.Checkbutton(root, text="Imaging Included (y/n)", variable=var, command = toggle_checkbox)
    checkbox.pack(side=TOP, anchor = "w", padx = 10)

    # FOLDER SELECTION
    folder_select_lbl = ttk.Label(root, text="Select folder to save images")
    # folder_select_lbl.pack(side=TOP, anchor = "w", padx = 10)

    folder_select_lbl_btn = ttk.Button(root, text="Choose folder", command=images_directory)
    # folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)

    dir_entry = ttk.Entry(root, width = 35)
    # dir_entry.pack(side=TOP, anchor = "w", padx = 10)


    # SAMPLE NAME

    sample_name_label = ttk.Label(root, text = "Sample Name")
    sample_name_label.pack(side=TOP, anchor = "w", padx = 10)
    sample_name_entry = ttk.Entry(root, width = 35)
    sample_name_entry.pack(side=TOP, anchor = "w", padx = 10)

    enter_sample_button = ttk.Button(root, text="Enter", command=read_sample_name)
    enter_sample_button.pack(side=TOP, anchor = "w", padx = 10, pady = 3)

    sample_prompt_text = ttk.Label(root, text = f"Enter cells for sample:")
    sample_prompt_text.pack(side=TOP, anchor = "w", padx = 10, pady = 3)

    # PLATE COUNTS
    plate_count_options = [1,2,3]
    plate_count_var = ttk.IntVar()
    plate_count_var.set(plate_count_options[0])
    plate_count_frame = ttk.Frame(root)
    plate_count_frame.pack(side=TOP, anchor = "w", padx=5)
    plate_count_label = ttk.Label(plate_count_frame, text = "Enter number of plates:")
    plate_count_label.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)
    plate_count_option_menu = ttk.OptionMenu(plate_count_frame, plate_count_var, *plate_count_options, command = plate_count_select)
    plate_count_option_menu.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)
    plate_count_enter = ttk.Button(plate_count_frame, text="Enter", command = disable_plate_count)
    plate_count_enter.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)

    plate_cells_frame = ttk.Frame(root)
    plate_cells_frame.pack(side=TOP, anchor = "w", padx=5)

    plate_cells_label = ttk.Label(plate_cells_frame, text = f"Enter cells for plate: ")
    plate_cells_label.pack(side=LEFT, anchor = "w", padx = 10, pady = 3)

    plate_cells_var = ttk.IntVar()
    plate_cells_var.set(1)
    plate_cells_entry = ttk.Entry(plate_cells_frame, textvariable=plate_cells_var)
    plate_cells_entry.pack(side=LEFT,anchor="w", padx=10,pady=5)

    # TABLE
    table = Table(root, rows=8, columns=12)
    table.pack(expand=True, fill = "both", side=TOP, anchor="w",padx = 10)

    enter_cells_frame = ttk.Frame(root)
    enter_cells_frame.pack(side=TOP, after=table, anchor = "w", padx=5)

    enter_cells = ttk.Button(enter_cells_frame, text="Enter", command=save_sample_cells)
    enter_cells.pack(side=LEFT, anchor = "w",pady=5)

    # SAVE TO JSON FILE

    save_to_json_btn = ttk.Button(root, text="Save to JSON", command=save_to_json)
    save_to_json_btn.pack(side=TOP, anchor = "s", padx = 10)

    # EXIT
    quit_btn = ttk.Button(root, text="Done", command=root.destroy)
    quit_btn.pack(side=TOP, padx = 10)

    root.mainloop()
