from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import tkfilebrowser
from table import *
import os
import subprocess
import numpy as np
import json
import platform

HOME_DIR = os.getcwd()
IMAGES_DIR = []
BULK_DIR = []
json_dict = {}
temp_json = {}
conditions_list = []
strains_list = []
condition_names = [""]
strain_names = [""]
all_conditions = []



# toggle_checkbox shows folder selection widgets if checkbox is checked
def toggle_checkbox_imaging():
    if var.get():
        folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        dir_entry.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        good_data_checkbox.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        dust_correction_checkbox.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        image_analysis_checkbox.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
    else:
        folder_select_lbl_btn.pack_forget()
        dir_entry.pack_forget()
        good_data_checkbox.pack_forget()
        dust_correction_checkbox.pack_forget()
        image_analysis_checkbox.pack_forget()

def send_directory():
    send_dir_entry.delete(0, END)
    send_dir_entry.insert(END, string=filedialog.askdirectory(initialdir="B:/"))

def images_directory():
    dir_entry.delete(0, "end")
    IMAGES_DIR.append(tkfilebrowser.askopendirnames(initialdir="C:/Users/Imaging Controller/Desktop/GEN5_IMAGE_LIBRARY"))
    for i in range(len(IMAGES_DIR)):
        if i == 0:
            dir_entry.insert(i, string=IMAGES_DIR[i][0]+",")
        else:
            dir_entry.insert(len(dir_entry.get())+1, string=IMAGES_DIR[i][0]+",")
    print(IMAGES_DIR)
    return IMAGES_DIR

def bulk_directory():
    bulk_dir_entry.delete(0, "end")
    BULK_DIR.append(filedialog.askopenfilenames(initialdir="C:/Users/Imaging Controller/Desktop/GEN5_IMAGE_LIBRARY"))
    for i in range(len(BULK_DIR)):
        if i == 0:
            bulk_dir_entry.insert(i, string=BULK_DIR[i][0]+",")
        else:
            bulk_dir_entry.insert(len(bulk_dir_entry.get())+1, string=BULK_DIR[i][0]+",")
    print(BULK_DIR)
    return BULK_DIR

def read_sample_and_strain_name():
    condition_names.clear()
    # sample_prompt_text.config(text = f"Enter cells for sample: {sample_name}")
    prompt_text.config(text = f"Enter cells for sample ({sample_name_entry.get()}) and plate ({str(plate_counter.get())})")
    enter_cells.configure(state="active")
    # plate_count_option_menu.configure(state="active")
    strain_names.clear()
    strain_name = strain_name_entry.get()

# save_sample_cells() saves the current list of sample_cells to
# the conditions_list dictionary, then clears the selected cells and the table
def save_sample_cells():
    nplates = plate_count_var.get()
    if len(conditions_list) < nplates:
        conditions_list.append({})
        strains_list.append({})
    curr_plate = plate_counter.get()
    sample_name = sample_name_entry.get()
    strain_name = strain_name_entry.get()
    conditions_list[curr_plate-1][sample_name] = sorted(list(sample_cells))
    strains_list[curr_plate-1][strain_name] = sorted(list(sample_cells))
    print(conditions_list)
    print(strains_list)
    if curr_plate < nplates:
        plate_counter.set(plate_counter.get() + 1)
    else:
        enter_cells.configure(state="disabled")
        plate_counter.set(1)
        sample_name = sample_name_entry.delete(0, ttk.END)
        strain_name = strain_name_entry.delete(0, ttk.END)
        # proceed_text = ttk.Label(text="Done! You can proceed to the next step.", fg="#70d000")
        # proceed_text.pack(side=LEFT, after=enter_cells, anchor="w",padx=5,pady=3)
    prompt_text.config(text = f"Enter cells for sample ({sample_name_entry.get()}) and plate ({str(plate_counter.get())})")
    sample_cells.clear()
    table.clear_color()
    

# disable_plate_count disables the plate_count_option_menu widget, prevents user from editing it after hitting "Enter"
def disable_plate_count():
    plate_count_option_menu.configure(state="disabled")
    
# def disable_media():
    # media_option_menu.configure(state="disabled")

def plate_count_select(value):
    plate_count_var.set(value)

def media_select(value):
    media_var.set(value)

# def save_plot_number():
#     plot_number = plot_number_entry.get("1.0", "end-1c")
#     with open("temp_plot_num.txt", "w") as fw:
#         fw.write(plot_number)

def create_new_window():
    json_dict["notes"] = notes_entry.get("1.0", "end-1c")
    json_dict["media"] = media_var.get()
    json_dict["acquisition_frequency"] = int(acquisition_freq.get("1.0", "end-1c"))
    print(IMAGES_DIR)
    json_dict["images_directory"] = [s[0].replace("\\", "/") for s in IMAGES_DIR]
    json_dict["bulk_data"] = [s[0] for s in BULK_DIR]
    json_dict["conditions"] = conditions_list
    json_dict["strains"] = strains_list
    if good_data_var.get():
        good_data_directory = "B:/Good imaging data"
    else:
        good_data_directory = ""
    if dust_var.get():
        dust_correction = "True"
    else:
        dust_correction = "False"
    if image_analysis_var.get():
        image_analysis = "True"
    else:
        image_analysis = "False"
    json_dict["dust_correction"] = dust_correction
    json_dict["good_data_directory"] = good_data_directory
    json_dict["image_analysis"] = image_analysis
    json_dict["experiment_directory"] = send_dir_entry.get()
    with open(HOME_DIR + "/temp_config.json", "w") as file:
        json.dump(json_dict, file, indent = 4)

    all_conditions = list(set(key for d in conditions_list for key in d))
    with open("temp_conditions.txt", "w") as fw:
        for cond in all_conditions:
            fw.write(cond)
            fw.write("\n")
    root.destroy()
    if platform.system() == "Windows":
        os.system('set LANG=en_US.UTF-8 && python3 ./GUI/plotoptions.py')
    elif platform.system() == "Darwin":
        os.system("export LANG='en_US.UTF-8' && python3 ./GUI/plotoptions.py")


def on_yes_click():
    file = open(HOME_DIR + "/temp_config.json", "r")
    temp_json = json.load(file)
    file.close()

    # repopulate notes entry
    notes_entry.insert(END, temp_json["notes"])
    
    # repopulate media option
    media_var.set(media_options[media_options.index(temp_json["media"])])
    plot_num_file = open(HOME_DIR + "/GUI/temp_plot_num.txt", "r")
    plot_number_entry.delete(1.0,END)
    plot_number_entry.insert(END, plot_num_file.read())
    plot_num_file.close()
    file.close()
    
    # acquisition frequency option
    acquisition_freq.delete(1.0, END)
    acquisition_freq.insert(END, temp_json["acquisition_frequency"])

    # experiment dir option
    send_dir_entry.delete(0, END)
    send_dir_entry.insert(END, temp_json["experiment_directory"])
    
    # image analysis checkbox option
    if temp_json["image_analysis"] == "True":
        
        var.set(1)
        toggle_checkbox_imaging()
        images_dirs = temp_json["images_directory"]
        for images_dir in images_dirs:
            IMAGES_DIR.append(tuple([images_dir]))
        # print([s[0].replace("\\", "/") for s in IMAGES_DIR])
        dir_entry.delete(0, END)
        dir_entry.insert(END, temp_json["images_directory"])
        if len(temp_json["good_data_directory"]) > 0:
            good_data_var.set(1)
        if temp_json["dust_correction"] == "True":
            dust_var.set(1)
        if temp_json["image_analysis"] == "True":
            image_analysis_var.set(1)
    bulk_dirs = temp_json["bulk_data"]
    for bulk_dir in bulk_dirs:
        BULK_DIR.append(tuple([bulk_dir]))
    bulk_dir_entry.delete(0,END)
    bulk_dir_entry.insert(END, temp_json["bulk_data"])
        
    
    
    
    
    root2.destroy()



if __name__ == "__main__":
    
    root = Tk()
    root.title("Scanner")
    root.geometry("1000x1000")
    
    
    
    frm = ttk.Frame(root)
    frm.grid()
    frm.pack(padx=2, pady=0)
    plate_counter = ttk.IntVar()
    plate_counter.set(1)
    # EXPERIMENT NOTES ENTRY
    notes_lbl = ttk.Label(root, text = "Experiment Notes")
    notes_lbl.pack(side=TOP, anchor = "w", padx = 10)

    notes_entry = Text(root, width = 35, height = 4)
    notes_entry.pack(side=TOP, anchor = "w", padx = 10)

    # EXPERIMENT METADATA
    media_options = ["LB", "LB+Glucose+CaCl2", "M9+Glucose+CA", "M9+Glucose+CA+2%Na",
                     "M9+Galactose+CA", "M9+Glycerol+CA", "M63+Glucose+CA"]
    media_var = ttk.StringVar()
    media_var.set(media_options[0])
    media_frame = ttk.Frame(root)
    media_frame.pack(side=TOP, anchor = "w", padx=5)
    media_label = ttk.Label(media_frame, text = "Enter media used for this experiment:")
    media_label.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)
    media_option_menu = ttk.OptionMenu(media_frame, media_var, *media_options, command = media_select)
    media_option_menu.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)
    # media_enter = ttk.Button(media_frame, text="Enter", command = disable_media)
    # media_enter.pack(side=LEFT, anchor = "w", padx = 5, pady = 3)

    # PLOT NUMBERS
    plot_number_frm = ttk.Frame(root)
    plot_number_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
    plot_number_lbl = ttk.Label(plot_number_frm, text = "Enter number of plots")
    plot_number_lbl.pack(side=LEFT, anchor = "w", padx = 10)
    plot_number_entry = Text(plot_number_frm, width = 5, height = 1)
    plot_number_entry.pack(side=LEFT, anchor = "w", padx = 10)
    # plot_number_btn = ttk.Button(plot_number_frm, text = "Enter", command = save_plot_number)
    # plot_number_btn.pack(side=LEFT, anchor = "w", padx = 10)

    # ACQUISITION FREQUENCY
    acquisition_lbl = ttk.Label(root, text = "Acquisition Frequency (#/hr)")
    acquisition_lbl.pack(side=TOP, anchor = "w", padx = 10)

    acquisition_freq = Text(root, width = 35, height = 1)
    acquisition_freq.pack(side=TOP, anchor = "w", padx = 10)

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

    # CHOOSE FOLDER TO SEND DATA
    send_folder_select_lbl_btn = ttk.Button(root, text="Choose folder to store this experiment", command=send_directory)
    send_dir_entry = ttk.Entry(root, width = 35)
    send_folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)
    send_dir_entry.pack(side=TOP, anchor = "w", padx = 10)

    # CHECKBOX FOR IMAGES DIRECTORY
    var = ttk.IntVar()
    checkbox = ttk.Checkbutton(root, text="Imaging included", variable=var, command = toggle_checkbox_imaging)
    checkbox.pack(side=TOP, anchor = "w", padx = 10)

    # FOLDER SELECTION
    images_dirs = []
    folder_select_lbl_btn = ttk.Button(root, text="Choose folder where images are located", command=images_directory)
    dir_entry = ttk.Entry(root, width = 35)

    # CHECKBOX FOR GOOD DATA
    good_data_var = ttk.IntVar()
    good_data_checkbox = ttk.Checkbutton(root, text="Good data", variable=good_data_var)

    # CHECKBOX FOR PERFORMING IMAGE ANALYSIS
    image_analysis_var = ttk.IntVar()
    image_analysis_checkbox = ttk.Checkbutton(root, text="Perform image analysis", variable=image_analysis_var)

    # CHECKBOX FOR DUST CORRECTION
    dust_var = ttk.IntVar()
    dust_correction_checkbox = ttk.Checkbutton(root, text="Perform dust correction", variable=dust_var)

    # BULK DATA SELECTION
    bulk_folder_select_lbl_btn = ttk.Button(root, text="Choose file for bulk data", command=bulk_directory)
    bulk_folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)
    bulk_dir_entry = ttk.Entry(root, width = 35)
    bulk_dir_entry.pack(side=TOP, anchor = "w", padx = 10)

    # SAMPLE NAME
    sample_name_label = ttk.Label(root, text = "Sample Name")
    sample_name_label.pack(side=TOP, anchor = "w", padx = 10)
    sample_name_entry = ttk.Entry(root, width = 35)
    sample_name_entry.pack(side=TOP, anchor = "w", padx = 10)
    # enter_sample_button = ttk.Button(root, text="Enter", command=read_sample_name)
    # enter_sample_button.pack(side=TOP, anchor = "w", padx = 10, pady = 3)

    # BRIDGES LAB STRAINS NAME
    strain_name_label = ttk.Label(root, text = "Lab Strain Name (with supplements)")
    strain_name_label.pack(side=TOP, anchor = "w", padx = 10)
    strain_name_entry = ttk.Entry(root, width = 35)
    strain_name_entry.pack(side=TOP, anchor = "w", padx = 10)
    enter_strain_button = ttk.Button(root, text="Enter", command=read_sample_and_strain_name)
    enter_strain_button.pack(side=TOP, anchor = "w", padx = 10, pady = 3)

    # ENTER CELLS
    # sample_prompt_text = ttk.Label(root, text = f"Enter cells for sample:")
    # sample_prompt_text.pack(side=TOP, anchor = "w", padx = 10, pady = 3)
    prompt_text = ttk.Label(root, text = f"Enter cells for sample () and plate ()")
    prompt_text.pack(side=TOP, anchor = "w", padx = 10, pady = 3)
    # plate_cells_frame = ttk.Frame(root)
    # plate_cells_frame.pack(side=TOP, anchor = "w", padx=5)
    # plate_cells_label = ttk.Label(plate_cells_frame, text = f"Enter cells for plate: ")
    # plate_cells_label.pack(side=LEFT, anchor = "w", padx = 10, pady = 3)
    # plate_cells_var = ttk.IntVar()
    # plate_cells_var.set(1)
    # plate_cells_entry = ttk.Entry(plate_cells_frame, textvariable=plate_cells_var)
    # plate_cells_entry.pack(side=LEFT,anchor="w", padx=10,pady=5)


    # TABLE
    table = Table(root, rows=9, columns=13)
    table.pack(expand=True, fill = "both", side=TOP, anchor="w",padx = 10)
    
    enter_cells_frame = ttk.Frame(root)
    enter_cells_frame.pack(side=LEFT, after=table, anchor = "w", padx=5)
    enter_cells = ttk.Button(enter_cells_frame, text="Enter", command=save_sample_cells)
    enter_cells.pack(side=LEFT, anchor = "w",pady=5)

    # NEXT
    next_btn = ttk.Button(root, text="Next", command=create_new_window)
    next_btn.pack(side=TOP, anchor = "e", padx = 10)
    
    # if not back_button_pressed.get():
    root2 = Tk()
    root2.title("Repopulate")
    root2.geometry("300x85")
    repopulate_frame = ttk.Frame(root2)
    repopulate_frame.pack(padx=2, pady=0)
    repopulate_label = ttk.Label(repopulate_frame, text="Repopulate values using config?")
    repopulate_label.pack(side=TOP, padx=4, pady=10)
    
    no_button = ttk.Button(repopulate_frame, text="No", fg="#bb2211",command=root2.destroy)
    no_button.pack(side=LEFT, after=repopulate_label, anchor = "sw", padx=4, pady=0)

    yes_button = ttk.Button(repopulate_frame, text="Yes",fg="#237911", command=on_yes_click)
    yes_button.pack(side=RIGHT, after=repopulate_label, anchor = "se", padx=4, pady=0)
    root2.mainloop()
    # else:
    #     root.destroy()
    
    root.mainloop()
