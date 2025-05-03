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
from natsort import natsorted
HOME_DIR = os.getcwd()
IMAGES_DIR = []
BULK_DIR = []
json_dict = {}
temp_json = {}

# toggle_checkbox shows folder selection widgets if checkbox is checked
def toggle_checkbox_imaging():
    if var.get():
        folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        dir_entry.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
        dust_correction_checkbox.pack(side=TOP, anchor = "w", padx = 10, before = bulk_folder_select_lbl_btn)
    else:
        folder_select_lbl_btn.pack_forget()
        dir_entry.pack_forget()
        dust_correction_checkbox.pack_forget()

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
    return IMAGES_DIR

def bulk_directory():
    bulk_dir_entry.delete(0, "end")
    BULK_DIR.append(filedialog.askopenfilenames(initialdir="C:/Users/Imaging Controller/Desktop/GEN5_IMAGE_LIBRARY"))
    for i in range(len(BULK_DIR)):
        if i == 0:
            bulk_dir_entry.insert(i, string=BULK_DIR[i][0]+",")
        else:
            bulk_dir_entry.insert(len(bulk_dir_entry.get())+1, string=BULK_DIR[i][0]+",")
    return BULK_DIR

def finish():
    json_dict["acquisition_frequency"] = int(acquisition_freq.get("1.0", "end-1c"))
    if upper_peak.get("1.0", "end-1c") == '':
        json_dict["upper_peak"] = 1.5
    else:
        json_dict["upper_peak"] = float(upper_peak.get("1.0", "end-1c"))
    if lower_peak.get("1.0", "end-1c") == '':
        json_dict["lower_peak"] = 0.5
    else:
        json_dict["lower_peak"] = float(lower_peak.get("1.0", "end-1c"))
    if upper_final.get("1.0", "end-1c") == '':
        json_dict["upper_final"] = 0.1
    else:
        json_dict["upper_final"] = float(upper_final.get("1.0", "end-1c"))
    json_dict["fixed_thresh"] = float(v1.get())
    json_dict["images_directory"] = [s[0].replace("\\", "/") for s in IMAGES_DIR]
    json_dict["bulk_data"] = [s[0] for s in BULK_DIR]
    if dust_var.get():
        dust_correction = "True"
    else:
        dust_correction = "False"
    json_dict["dust_correction"] = dust_correction
    json_dict["experiment_directory"] = send_dir_entry.get()
    with open(HOME_DIR + "/experiment_config.json", "w") as file:
        json.dump(json_dict, file, indent = 4)
    root.destroy()

if __name__ == "__main__":

    root = Tk()
    root.title("Scanner")
    root.geometry("1000x500")
    frm = ttk.Frame(root)
    frm.grid()
    frm.pack(padx=2, pady=0)
    plate_counter = ttk.IntVar()
    plate_counter.set(1)

    # ACQUISITION FREQUENCY
    acquisition_lbl = ttk.Label(root, text = "Acquisition Frequency (#/hr)")
    acquisition_lbl.pack(side=TOP, anchor = "w", padx = 10)

    acquisition_freq = Text(root, width = 35, height = 1)
    acquisition_freq.pack(side=TOP, anchor = "w", padx = 10)

    # SLIDER FOR SETTING THRESHOLD
    v1 = DoubleVar()
    s1 = Scale(root, variable = v1,
               from_ = 0.0000, to = 1.0000, resolution=0.0001,
               orient = HORIZONTAL)
    v1.set(0.04)
    l3 = Label(root, text = "Threshold")
    s1.pack(anchor = "w")
    l3.pack(anchor = "w")

    # UPPER, LOWER FOR PEAK
    upper_peak_lbl = ttk.Label(root, text = "Upper FC for peak relative to mean peak for plate (default = 1.5)")
    upper_peak_lbl.pack(side=TOP, anchor = "w", padx = 10)
    upper_peak = Text(root, width = 35, height = 1)
    upper_peak.pack(side=TOP, anchor = "w", padx = 10)

    lower_peak_lbl = ttk.Label(root, text = "Lower FC for peak relative to mean peak for plate (default = 0.5)")
    lower_peak_lbl.pack(side=TOP, anchor = "w", padx = 10)
    lower_peak = Text(root, width = 35, height = 1)
    lower_peak.pack(side=TOP, anchor = "w", padx = 10)

    # UPPER FOR FINAL TIMEPOINT
    upper_final_lbl = ttk.Label(root, text = "Upper FC of final biomass relative to mean peak of the plate (default = 0.1)")
    upper_final_lbl.pack(side=TOP, anchor = "w", padx = 10)
    upper_final = Text(root, width = 35, height = 1)
    upper_final.pack(side=TOP, anchor = "w", padx = 10)

    # CHOOSE FOLDER TO SEND DATA
    send_folder_select_lbl_btn = ttk.Button(root, text="Choose folder to store this experiment", command=send_directory)
    send_dir_entry = ttk.Entry(root, width = 35)
    send_folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)
    send_dir_entry.pack(side=TOP, anchor = "w", padx = 10)

    # FOLDER SELECTION
    images_dirs = []
    folder_select_lbl_btn = ttk.Button(root, text="Choose folder where images are located", command=images_directory)
    folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)
    dir_entry = ttk.Entry(root, width = 35)
    dir_entry.pack(side=TOP, anchor = "w", padx = 10)

    # CHECKBOX FOR DUST CORRECTION
    dust_var = ttk.IntVar()
    dust_correction_checkbox = ttk.Checkbutton(root, text="Perform dust correction", variable=dust_var)

    # BULK DATA SELECTION
    bulk_folder_select_lbl_btn = ttk.Button(root, text="Choose file for bulk data", command=bulk_directory)
    bulk_folder_select_lbl_btn.pack(side=TOP, anchor = "w", padx = 10)
    bulk_dir_entry = ttk.Entry(root, width = 35)
    bulk_dir_entry.pack(side=TOP, anchor = "w", padx = 10)

    # NEXT
    next_btn = ttk.Button(root, text="Next", command=finish)
    next_btn.pack(side=TOP, anchor = "e", padx = 10)

    root.mainloop()
