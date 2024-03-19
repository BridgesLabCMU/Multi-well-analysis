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

plot_types = {"plot1": "", "plot2": "", "plot3":""}
plot_dtypes = {"plot1":[], "plot2":[], "plot3":[]}
plot_xaxes = {"plot1":"", "plot2":"", "plot3":""}
plot_normalizations = {"plot1":"", "plot2":"", "plot3":""}
plot_titles = {"plot1":"", "plot2":"", "plot3":""}
plot_ylabs = {"plot1":"", "plot2":"", "plot3":""}
plot_xlabs = {"plot1":"", "plot2":"", "plot3":""}
plot_filenames = {"plot1":"", "plot2":"", "plot3":""}
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
    plot1_normalization_option_menu['menu'].add_command(label=sample_name, command=ttk._setit(plot1_normalization_selection, sample_name))
    plot2_normalization_option_menu['menu'].add_command(label=sample_name, command=ttk._setit(plot2_normalization_selection, sample_name))
    plot3_normalization_option_menu['menu'].add_command(label=sample_name, command=ttk._setit(plot3_normalization_selection, sample_name))

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
    store_plot_options()
    json_dict["plot_types"] = plot_types
    json_dict["plot_dtypes"] = plot_dtypes
    json_dict["plot_xaxis"] = plot_xaxes
    json_dict["plot_normalization"] = plot_normalizations
    json_dict["plot_titles"] = plot_titles
    json_dict["plot_xlabs"] = plot_xlabs
    json_dict["plot_ylabs"] = plot_ylabs
    json_dict["plot_filenames"] = plot_filenames
    json_dict["sig"] = 2
    json_dict["blockDiameter"] = [501, 101]
    json_dict["shift_thresh"] = 50
    with open(HOME_DIR + "/test_configs.json", "w") as file:
        json.dump(json_dict, file, indent = 4)

# calling plot option storing functions
def store_plot_options():
    save_plot_types(plot1_type_selection, plot2_type_selection, plot3_type_selection)
    save_plot_dtypes(plot1_dtype_selection, plot2_dtype_selection, plot3_dtype_selection)
    save_plot_xaxis(plot1_xaxis_selection, plot2_xaxis_selection, plot3_xaxis_selection)
    save_plot_normalization(plot1_normalization_selection, plot2_normalization_selection, plot3_normalization_selection)
    save_plot_titles(plot1_title_entry, plot2_title_entry, plot3_title_entry)
    save_plot_xlabs(plot1_xlabs_entry, plot2_xlabs_entry, plot3_xlabs_entry)
    save_plot_ylabs(plot1_ylabs_entry, plot2_ylabs_entry, plot3_ylabs_entry)
    save_plot_filenames(plot1_filename_entry, plot2_filename_entry, plot3_filename_entry)

# storing plot options in list variables
def save_plot_types(plot1_type_selection, plot2_type_selection, plot3_type_selection):
    plot_types["plot1"] = plot1_type_selection.get()
    plot_types["plot2"] = plot2_type_selection.get()
    plot_types["plot3"] = plot3_type_selection.get()
def save_plot_dtypes(plot1_dtype_selection, plot2_dtype_selection, plot3_dtype_selection):
    plot_dtypes["plot1"] = plot1_dtype_selection.get().split(",")
    plot_dtypes["plot2"] = plot2_dtype_selection.get().split(",")
    plot_dtypes["plot3"] = plot3_dtype_selection.get().split(",")
def save_plot_xaxis(plot1_xaxis_selection, plot2_xaxis_selection, plot3_xaxis_selection):
    plot_xaxes["plot1"] = plot1_xaxis_selection.get()
    plot_xaxes["plot2"] = plot2_xaxis_selection.get()
    plot_xaxes["plot3"] = plot3_xaxis_selection.get()
def save_plot_normalization(plot1_normalization_selection, plot2_normalization_selection, plot3_normalization_selection):
    plot_normalizations["plot1"] = plot1_normalization_selection.get()
    plot_normalizations["plot2"] = plot2_normalization_selection.get()
    plot_normalizations["plot3"] = plot3_normalization_selection.get()
def save_plot_titles(plot1_title_entry, plot2_title_entry, plot3_title_entry):
    plot_titles["plot1"] = plot1_title_entry.get()
    plot_titles["plot2"] = plot2_title_entry.get()
    plot_titles["plot3"] = plot3_title_entry.get()
def save_plot_xlabs(plot1_xlabs_entry, plot2_xlabs_entry, plot3_xlabs_entry):
    plot_xlabs["plot1"] = plot1_xlabs_entry.get()
    plot_xlabs["plot2"] = plot2_xlabs_entry.get()
    plot_xlabs["plot3"] = plot3_xlabs_entry.get()
def save_plot_ylabs(plot1_ylabs_entry, plot2_ylabs_entry, plot3_ylabs_entry):
    plot_ylabs["plot1"] = plot1_ylabs_entry.get()
    plot_ylabs["plot2"] = plot2_ylabs_entry.get()
    plot_ylabs["plot3"] = plot3_ylabs_entry.get()
def save_plot_filenames(plot1_filename_entry, plot2_filename_entry, plot3_filename_entry):
    plot_filenames["plot1"] = plot1_filename_entry.get()
    plot_filenames["plot2"] = plot2_filename_entry.get()
    plot_filenames["plot3"] = plot3_filename_entry.get()
def plate_count_select(value):
    plate_count_var.set(value)
def plot1_type_select(value):
    plot1_type_selection.set(value)
def plot2_type_select(value):
    plot2_type_selection.set(value)
def plot3_type_select(value):
    plot3_type_selection.set(value)
def plot1_dtype_select(value):
    plot1_dtype_selection.set(value)
def plot2_dtype_select(value):
    plot2_dtype_selection.set(value)
def plot3_dtype_select(value):
    plot3_dtype_selection.set(value)
def plot1_xaxis_select(value):
    plot1_xaxis_selection.set(value)
def plot2_xaxis_select(value):
    plot2_xaxis_selection.set(value)
def plot3_xaxis_select(value):
    plot3_xaxis_selection.set(value)
def plot1_normalization_select(value):
    plot1_normalization_selection.set(value)
def plot2_normalization_select(value):
    plot2_normalization_selection.set(value)
def plot3_normalization_select(value):
    plot3_normalization_selection.set(value)

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

    # PLOT OPTIONS DROP DOWN MENUS

    plot_type_options_frame = ttk.Frame(root)
    plot_type_options_frame.pack(side=TOP, after=enter_cells_frame, anchor = "w", padx=5)
    plot_type_options = ["line", "jitter", "heatmap", ""]
    plot1_type_selection = ttk.StringVar()
    plot2_type_selection = ttk.StringVar()
    plot3_type_selection = ttk.StringVar()

    plot1_type_selection.set(plot_type_options[0])
    plot2_type_selection.set(plot_type_options[1])
    plot3_type_selection.set(plot_type_options[2])

    plot_type_selection_label = ttk.Label(plot_type_options_frame, text = "Select plot types for plot1, plot2, and plot3:")
    plot_type_selection_label.pack(side=LEFT, anchor = "nw", padx=5, pady=5)

    plot1_type_option_menu = ttk.OptionMenu(plot_type_options_frame, plot1_type_selection, *plot_type_options, command = plot1_type_select)
    plot1_type_option_menu.pack(side=LEFT, after = plot_type_selection_label, anchor = "w", padx = 5)

    plot2_type_option_menu = ttk.OptionMenu(plot_type_options_frame, plot2_type_selection, *plot_type_options, command = plot2_type_select)
    plot2_type_option_menu.pack(side=LEFT, after = plot1_type_option_menu, anchor = "w", padx = 5)

    plot3_type_option_menu = ttk.OptionMenu(plot_type_options_frame, plot3_type_selection, *plot_type_options, command = plot3_type_select)
    plot3_type_option_menu.pack(side=LEFT, after = plot2_type_option_menu, anchor = "w", padx = 5)



    # PLOT DTYPE OPTIONS DROP DOWN MENUS
    plot_dtype_options_frame = ttk.Frame(root)
    plot_dtype_options_frame.pack(side=TOP, anchor = "w", padx=5, pady=5)

    plot_dtype_options = ["lum, OD", "BF_imaging", "CY5_imaging"]
    plot1_dtype_selection = ttk.StringVar()
    plot2_dtype_selection = ttk.StringVar()
    plot3_dtype_selection = ttk.StringVar()
    plot1_dtype_selection.set(plot_dtype_options[0])
    plot2_dtype_selection.set(plot_dtype_options[1])
    plot3_dtype_selection.set(plot_dtype_options[2])


    plot_dtype_selection_label = ttk.Label(plot_dtype_options_frame, text = "Select plot datatypes for plot1, plot2, and plot3:")
    plot_dtype_selection_label.pack(side=LEFT, anchor = "nw", padx=5, pady=5)

    plot1_dtype_option_menu = ttk.OptionMenu(plot_dtype_options_frame, plot1_dtype_selection, *plot_dtype_options, command = plot1_dtype_select)
    plot1_dtype_option_menu.pack(side=LEFT, after = plot_dtype_selection_label, anchor = "w", padx = 5)

    plot2_dtype_option_menu = ttk.OptionMenu(plot_dtype_options_frame, plot2_dtype_selection, *plot_dtype_options, command = plot2_dtype_select)
    plot2_dtype_option_menu.pack(side=LEFT, after = plot1_dtype_option_menu, anchor = "w", padx = 5)

    plot3_dtype_option_menu = ttk.OptionMenu(plot_dtype_options_frame, plot3_dtype_selection, *plot_dtype_options, command = plot3_dtype_select)
    plot3_dtype_option_menu.pack(side=LEFT, after = plot2_dtype_option_menu, anchor = "w", padx = 5)

    # PLOT xaxis OPTIONS
    plot1_xaxis_selection = ttk.StringVar()
    plot2_xaxis_selection = ttk.StringVar()
    plot3_xaxis_selection = ttk.StringVar()
    plot1_xaxis_selection.set("")
    plot2_xaxis_selection.set("")
    plot3_xaxis_selection.set("")

    plot_xaxis_frm = ttk.Frame(root)
    plot_xaxis_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)

    plot1_xaxis_label = ttk.Label(plot_xaxis_frm, text="Plot 1 xaxis:")
    plot1_xaxis_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_xaxis_option_menu = ttk.OptionMenu(plot_xaxis_frm, plot1_xaxis_selection, *condition_names, command = plot1_xaxis_select)
    plot1_xaxis_option_menu.pack(side=LEFT, after = plot1_xaxis_label, anchor = "w", padx = 5)

    plot2_xaxis_label = ttk.Label(plot_xaxis_frm, text="Plot 2 xaxis:")
    plot2_xaxis_label.pack(side=LEFT, after = plot1_xaxis_option_menu, anchor = "w", padx = 5)
    plot2_xaxis_option_menu = ttk.OptionMenu(plot_xaxis_frm, plot2_xaxis_selection, *condition_names, command = plot2_xaxis_select)
    plot2_xaxis_option_menu.pack(side=LEFT, after = plot2_xaxis_label, anchor = "w", padx = 5)

    plot3_xaxis_label = ttk.Label(plot_xaxis_frm, text="Plot 3 xaxis:")
    plot3_xaxis_label.pack(side=LEFT, after = plot2_xaxis_option_menu, anchor = "w", padx = 5)
    plot3_xaxis_option_menu = ttk.OptionMenu(plot_xaxis_frm, plot3_xaxis_selection, *condition_names, command = plot3_xaxis_select)
    plot3_xaxis_option_menu.pack(side=LEFT, after = plot3_xaxis_label, anchor = "w", padx = 5)


    # PLOT NORMALIZATION OPTIONS
    plot1_normalization_selection = ttk.StringVar()
    plot2_normalization_selection = ttk.StringVar()
    plot3_normalization_selection = ttk.StringVar()
    plot1_normalization_selection.set("")
    plot2_normalization_selection.set("")
    plot3_normalization_selection.set("")

    plot_normalization_frm = ttk.Frame(root)
    plot_normalization_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)

    plot1_normalization_label = ttk.Label(plot_normalization_frm, text="Plot 1 normalization:")
    plot1_normalization_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_normalization_option_menu = ttk.OptionMenu(plot_normalization_frm, plot1_normalization_selection, *condition_names, command = plot1_normalization_select)
    plot1_normalization_option_menu.pack(side=LEFT, after = plot1_normalization_label, anchor = "w", padx = 5)

    plot2_normalization_label = ttk.Label(plot_normalization_frm, text="Plot 2 normalization:")
    plot2_normalization_label.pack(side=LEFT, after = plot1_normalization_option_menu, anchor = "w", padx = 5)
    plot2_normalization_option_menu = ttk.OptionMenu(plot_normalization_frm, plot2_normalization_selection, *condition_names, command = plot2_normalization_select)
    plot2_normalization_option_menu.pack(side=LEFT, after = plot2_normalization_label, anchor = "w", padx = 5)

    plot3_normalization_label = ttk.Label(plot_normalization_frm, text="Plot 3 normalization:")
    plot3_normalization_label.pack(side=LEFT, after = plot2_normalization_option_menu, anchor = "w", padx = 5)
    plot3_normalization_option_menu = ttk.OptionMenu(plot_normalization_frm, plot3_normalization_selection, *condition_names, command = plot3_normalization_select)
    plot3_normalization_option_menu.pack(side=LEFT, after = plot3_normalization_label, anchor = "w", padx = 5)


    # PLOT TITLES
    plot_title_frm = ttk.Frame(root)
    plot_title_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)


    plot1_title_label = ttk.Label(plot_title_frm, text="Plot 1 title:")
    plot1_title_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_title_entry = ttk.Entry(plot_title_frm)
    plot1_title_entry.pack(side=LEFT, after = plot1_title_label, anchor = "w", padx = 5)

    plot2_title_label = ttk.Label(plot_title_frm, text="Plot 2 title:")
    plot2_title_label.pack(side=LEFT, after = plot1_title_entry, anchor = "w", padx = 5)
    plot2_title_entry = ttk.Entry(plot_title_frm)
    plot2_title_entry.pack(side=LEFT, after = plot2_title_label, anchor = "w", padx = 5)

    plot3_title_label = ttk.Label(plot_title_frm, text="Plot 3 title:")
    plot3_title_label.pack(side=LEFT, after = plot2_title_entry, anchor = "w", padx = 5)
    plot3_title_entry = ttk.Entry(plot_title_frm)
    plot3_title_entry.pack(side=LEFT, after = plot3_title_label, anchor = "w", padx = 5)

    # PLOT YLABS
    plot_ylabs_frm = ttk.Frame(root)
    plot_ylabs_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)


    plot1_ylabs_label = ttk.Label(plot_ylabs_frm, text="Plot 1 y-label:")
    plot1_ylabs_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_ylabs_entry = ttk.Entry(plot_ylabs_frm)
    plot1_ylabs_entry.pack(side=LEFT, after = plot1_ylabs_label, anchor = "w", padx = 5)

    plot2_ylabs_label = ttk.Label(plot_ylabs_frm, text="Plot 2 y-label:")
    plot2_ylabs_label.pack(side=LEFT, after = plot1_ylabs_entry, anchor = "w", padx = 5)
    plot2_ylabs_entry = ttk.Entry(plot_ylabs_frm)
    plot2_ylabs_entry.pack(side=LEFT, after = plot2_ylabs_label, anchor = "w", padx = 5)

    plot3_ylabs_label = ttk.Label(plot_ylabs_frm, text="Plot 3 y-label:")
    plot3_ylabs_label.pack(side=LEFT, after = plot2_ylabs_entry, anchor = "w", padx = 5)
    plot3_ylabs_entry = ttk.Entry(plot_ylabs_frm)
    plot3_ylabs_entry.pack(side=LEFT, after = plot3_ylabs_label, anchor = "w", padx = 5)



    # PLOT XLABS
    plot_xlabs_frm = ttk.Frame(root)
    plot_xlabs_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)


    plot1_xlabs_label = ttk.Label(plot_xlabs_frm, text="Plot 1 x-label:")
    plot1_xlabs_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_xlabs_entry = ttk.Entry(plot_xlabs_frm)
    plot1_xlabs_entry.pack(side=LEFT, after = plot1_xlabs_label, anchor = "w", padx = 5)

    plot2_xlabs_label = ttk.Label(plot_xlabs_frm, text="Plot 2 x-label:")
    plot2_xlabs_label.pack(side=LEFT, after = plot1_xlabs_entry, anchor = "w", padx = 5)
    plot2_xlabs_entry = ttk.Entry(plot_xlabs_frm)
    plot2_xlabs_entry.pack(side=LEFT, after = plot2_xlabs_label, anchor = "w", padx = 5)

    plot3_xlabs_label = ttk.Label(plot_xlabs_frm, text="Plot 3 x-label:")
    plot3_xlabs_label.pack(side=LEFT, after = plot2_xlabs_entry, anchor = "w", padx = 5)
    plot3_xlabs_entry = ttk.Entry(plot_xlabs_frm)
    plot3_xlabs_entry.pack(side=LEFT, after = plot3_xlabs_label, anchor = "w", padx = 5)


    # PLOT FILENAMES
    plot_filename_frm = ttk.Frame(root)
    plot_filename_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)


    plot1_filename_label = ttk.Label(plot_filename_frm, text="Plot 1 filename:")
    plot1_filename_label.pack(side=LEFT, anchor = "w", padx = 5)
    plot1_filename_entry = ttk.Entry(plot_filename_frm)
    plot1_filename_entry.pack(side=LEFT, after = plot1_filename_label, anchor = "w", padx = 5)

    plot2_filename_label = ttk.Label(plot_filename_frm, text="Plot 2 filename:")
    plot2_filename_label.pack(side=LEFT, after = plot1_filename_entry, anchor = "w", padx = 5)
    plot2_filename_entry = ttk.Entry(plot_filename_frm)
    plot2_filename_entry.pack(side=LEFT, after = plot2_filename_label, anchor = "w", padx = 5)

    plot3_filename_label = ttk.Label(plot_filename_frm, text="Plot 3 filename:")
    plot3_filename_label.pack(side=LEFT, after = plot2_filename_entry, anchor = "w", padx = 5)
    plot3_filename_entry = ttk.Entry(plot_filename_frm)
    plot3_filename_entry.pack(side=LEFT, after = plot3_filename_label, anchor = "w", padx = 5)



    # SAVE TO JSON FILE

    save_to_json_btn = ttk.Button(root, text="Save to JSON", command=save_to_json)
    save_to_json_btn.pack(side=TOP, anchor = "s", padx = 10)



    # EXIT
    quit_btn = ttk.Button(root, text="Done", command=root.destroy)
    quit_btn.pack(side=TOP, padx = 10)

    # empty_label = ttk.Label(root, text="")
    # empty_label.pack(side = TOP, padx = 10)

    root.mainloop()
