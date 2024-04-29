from tkinter import ttk
import tkinter as tk
from tkinter import *
from tkinter.constants import *
from gui import *
import json
import numpy as np
from PIL import Image, ImageTk

class VerticalScrolledFrame(ttk.Frame):
    def __init__(self, parent, *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)

        # Create a canvas object and a vertical scrollbar for scrolling it.
        vscrollbar = ttk.Scrollbar(self, orient=VERTICAL)
        vscrollbar.pack(fill=Y, side=RIGHT, expand=FALSE)
        canvas = tk.Canvas(self, width=2000, height=1000, bd=0, highlightthickness=0,
                           yscrollcommand=vscrollbar.set)
        canvas.pack(side=LEFT, fill=BOTH, expand=TRUE)
        vscrollbar.config(command=canvas.yview)

        # Reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # Create a frame inside the canvas which will be scrolled with it.
        self.interior = interior = ttk.Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor=NW)

        # Track changes to the canvas and frame width and sync them,
        # also updating the scrollbar.
        def _configure_interior(event):
            # Update the scrollbars to match the size of the inner frame.
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # Update the canvas's width to fit the inner frame.
                canvas.config(width=interior.winfo_reqwidth())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # Update the inner frame's width to fill the canvas.
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)

num_plots = 0
with open("temp_plot_num.txt", "r") as fr:
    num_plots = int(fr.read())
    fr.close()

condition_names = []
with open("temp_conditions.txt", "r") as fr:
    for line in fr:
        if line[:-1] not in condition_names:
            condition_names.append(line[:-1])
    fr.close()

plot_types = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_dtypes = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_numerators = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_denominators = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_conditions = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_normalizations = {f"plot{i}": "" for i in range(1, num_plots+1)}
normalization_methods = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_xaxes = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_scales = {f"plot{i}": [] for i in range(1, num_plots+1)}
plot_colors = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_fonts = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_titles = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_ylabs = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_xlabs = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_xticks = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_yticks = {f"plot{i}": "" for i in range(1, num_plots+1)}
color_labels = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_sizes = {f"plot{i}": "" for i in range(1, num_plots+1)}
ylims = {f"plot{i}": "" for i in range(1, num_plots+1)}
dose_concs = {f"plot{i}": "" for i in range(1, num_plots+1)}
plot_filenames = {f"plot{i}": "" for i in range(1, num_plots+1)}


if __name__ == "__main__":
    class PlotOptions(tk.Tk):
        def __init__(self, *args, **kwargs):
            def save_to_json():
                with open("temp_config.json", "r") as file:
                    json_dict = json.load(file)
                store_plot_options()
                json_dict["plot_types"] = plot_types
                json_dict["plot_dtypes"] = plot_dtypes
                json_dict["plot_numerators"] = plot_numerators
                json_dict["plot_denominators"] = plot_denominators
                json_dict["plot_conditions"] = plot_conditions
                json_dict["plot_normalization"] = plot_normalizations
                json_dict["normalization_method"] = normalization_methods
                json_dict["plot_xaxis"] = plot_xaxes
                json_dict["plot_scale"] = plot_scales
                json_dict["plot_color"] = plot_colors
                json_dict["plot_font"] = plot_fonts
                json_dict["plot_titles"] = plot_titles
                json_dict["plot_xlabs"] = plot_xlabs
                json_dict["plot_ylabs"] = plot_ylabs
                json_dict["plot_xticks"] = plot_xticks
                json_dict["plot_yticks"] = plot_yticks
                json_dict["color_label"] = color_labels
                json_dict["ylims"] = ylims
                json_dict["plot_size"] = plot_sizes
                json_dict["dose_concs"] = dose_concs
                json_dict["plot_filenames"] = plot_filenames
                json_dict["sig"] = 2
                json_dict["blockDiameter"] = [501, 101]
                json_dict["shift_thresh"] = 50
                with open(HOME_DIR + "/experiment_config.json", "w") as file:
                    json.dump(json_dict, file, indent = 4)
            def store_plot_options():
                save_plot_types(plot_type_selections)
                save_plot_dtypes(plot_dtype_selections)
                save_plot_numerators(plot_numerator_selections)
                save_plot_denominators(plot_denominator_selections)
                save_plot_conditions(plot_condition_selections)
                save_plot_normalization(plot_normalization_selections)
                save_normalization_methods(normalization_method_selections)
                save_plot_xaxis(plot_xaxis_selections)
                save_plot_scale(plot_scale_selections)
                save_plot_color(plot_color_entries)
                save_plot_font(plot_font_selections)
                save_plot_titles(plot_title_entries)
                save_plot_xlabs(plot_xlab_entries)
                save_plot_ylabs(plot_ylab_entries)
                save_plot_xticks(plot_xtick_entries)
                save_plot_yticks(plot_ytick_entries)
                save_color_labels(color_label_entries)
                save_plot_sizes(plot_size_entries)
                save_ylims(ylim_entries)
                save_dose_concs(dose_conc_entries)
                save_plot_filenames(plot_filename_entries)

            def save_plot_types(plot_type_selections):
                for i in range(0,len(plot_type_selections)):
                    plot_types[f"plot{i+1}"] = plot_type_selections[i].get()
            def save_plot_dtypes(plot_dtype_selections):
                for i in range(0,len(plot_dtype_selections)):
                    plot_dtypes[f"plot{i+1}"] = plot_dtype_selections[i]
            def save_plot_numerators(plot_numerator_selections):
                for i in range(0,len(plot_numerator_selections)):
                    plot_numerators[f"plot{i+1}"] = plot_numerator_selections[i]
            def save_plot_denominators(plot_denominator_selections):
                for i in range(0,len(plot_denominator_selections)):
                    plot_denominators[f"plot{i+1}"] = plot_denominator_selections[i]
            def save_plot_conditions(plot_condition_selections):
                for i in range(0, len(plot_condition_selections)):
                    plot_conditions[f"plot{i+1}"] = plot_condition_selections[i]
            def save_plot_normalization(plot_normalization_selections):
                for i in range(0,len(plot_normalization_selections)):
                    plot_normalizations[f"plot{i+1}"] = plot_normalization_selections[i].get()
            def save_normalization_methods(normalization_method_selections):
                for i in range(0,len(normalization_method_selections)):
                    normalization_methods[f"plot{i+1}"] = normalization_method_selections[i]
            def save_plot_xaxis(plot_xaxis_selections):
                for i in range(0,len(plot_xaxis_selections)):
                    plot_xaxes[f"plot{i+1}"] = plot_xaxis_selections[i].get()
            def save_plot_scale(plot_scale_selections):
                for i in range(0,len(plot_scale_selections)):
                    plot_scales[f"plot{i+1}"] = plot_scale_selections[i]
            def save_plot_color(plot_color_entries):
                for i in range(0, len(plot_color_entries)):
                    plot_colors[f"plot{i+1}"] = plot_color_entries[i].get()
            def save_plot_font(plot_font_selections):
                for i in range(0,len(plot_font_selections)):
                    plot_fonts[f"plot{i+1}"] = plot_font_selections[i].get()
            def save_plot_titles(plot_title_entries):
                for i in range(0, len(plot_title_entries)):
                    plot_titles[f"plot{i+1}"] = plot_title_entries[i].get()
            def save_plot_xlabs(plot_xlab_entries):
                for i in range(0, len(plot_xlab_entries)):
                    plot_xlabs[f"plot{i+1}"] = plot_xlab_entries[i].get()
            def save_plot_ylabs(plot_ylab_entries):
                for i in range(0, len(plot_ylab_entries)):
                    plot_ylabs_i = plot_ylab_entries[i].get()
                    if plot_ylabs_i != "":
                        if "," in plot_ylabs_i:
                            plot_ylabs[f"plot{i+1}"] = plot_ylabs_i.split(",")
                        else:
                            plot_ylabs[f"plot{i+1}"] = [plot_ylabs_i]
                    else:
                        plot_ylabs[f"plot{i+1}"] = [plot_ylabs_i]
            def save_plot_xticks(plot_xtick_entries):
                for i in range(0, len(plot_xtick_entries)):
                    plot_xticks_i = plot_xtick_entries[i].get()
                    if plot_xticks_i != "":
                        if "," in plot_xticks_i:
                            plot_xticks[f"plot{i+1}"] = plot_xticks_i.split(",")
                    else:
                        plot_xticks[f"plot{i+1}"] = [plot_xticks_i]
            def save_plot_yticks(plot_ytick_entries):
                for i in range(0, len(plot_ytick_entries)):
                    plot_yticks_i = plot_ytick_entries[i].get()
                    if plot_yticks_i != "":
                        if "," in plot_yticks_i:
                            plot_yticks[f"plot{i+1}"] = plot_yticks_i.split(",")
                    else:
                        plot_yticks[f"plot{i+1}"] = [plot_yticks_i]
            def save_color_labels(color_label_entries):
                for i in range(0, len(color_label_entries)):
                    color_labels[f"plot{i+1}"] = color_label_entries[i].get()
            def save_plot_sizes(plot_size_entries):
                for i in range(0, len(plot_size_entries)):
                    plot_sizes_i = plot_size_entries[i].get()
                    if plot_sizes_i != "":
                        if "," in plot_sizes_i:
                            plot_sizes[f"plot{i+1}"] = [int(x) for x in plot_sizes_i.split(",")]
                    else:
                        plot_sizes[f"plot{i+1}"] = [300, 250]
            def save_ylims(ylim_entries):
                for i in range(0, len(ylim_entries)):
                    ylims_i = ylim_entries[i].get()
                    if ylims_i != "":
                        if "," in ylims_i:
                            ylim_list = []
                            for x in ylims_i.split(","):
                                if x == "":
                                    ylim_list.append("nothing")
                                else:
                                    ylim_list.append(int(x))
                            ylims[f"plot{i+1}"] = ylim_list
                    else:
                        ylims[f"plot{i+1}"] = "default"
            def save_dose_concs(dose_conc_entries):
                for i in range(0, len(dose_conc_entries)):
                    dose_concs_i = dose_conc_entries[i].get()
                    if dose_concs_i != "":
                        if "," in dose_concs_i:
                            dose_concs[f"plot{i+1}"] = [int(x) for x in dose_concs_i.split(",")]
                    else:
                        dose_concs[f"plot{i+1}"] = []
            def save_plot_filenames(plot_filename_entries):
                for i in range(0, len(plot_filename_entries)):
                    plot_filenames[f"plot{i+1}"] = plot_filename_entries[i].get()

            def plot_type_select(value, i):
                plot_type_selections[i-1].set(value)
            def plot_normalization_select(value, i):
                plot_normalization_selections[i-1].set(value)
            def plot_xaxis_select(value, i):
                plot_xaxis_selections[i-1].set(value)
            def plot_font_select(value, i):
                plot_font_selections[i-1].set(value)

            def normalization_method_listbox_on_select(event):
                normalization_method_selections.clear()
                for i, listbox in enumerate(normalization_method_listboxes):
                    selected_indices = listbox.curselection()
                    if len(normalization_method_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in normalization_method_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: normalization_method_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in normalization_method_prev_selected[i]:
                                items.append(listbox.get(idx))
                        normalization_method_selections.append(items)
                    else:
                        normalization_method_selections.append([listbox.get(idx) for idx in selected_indices])
                normalization_method_prev_selected.clear()
                normalization_method_prev_selected.extend(normalization_method_selections)

            def plot_scale_listbox_on_select(event):
                plot_scale_selections.clear()
                for i, listbox in enumerate(plot_scale_listboxes):
                    selected_indices = listbox.curselection()
                    if len(plot_scale_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in plot_scale_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: plot_scale_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in plot_scale_prev_selected[i]:
                                items.append(listbox.get(idx))
                        plot_scale_selections.append(items)
                    else:
                        plot_scale_selections.append([listbox.get(idx) for idx in selected_indices])
                plot_scale_prev_selected.clear()
                plot_scale_prev_selected.extend(plot_scale_selections)

            def dtype_listbox_on_select(event):
                plot_dtype_selections.clear()
                for i, listbox in enumerate(plot_dtype_listboxes):
                    selected_indices = listbox.curselection()
                    if len(plot_dtype_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in plot_dtype_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: plot_dtype_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in plot_dtype_prev_selected[i]:
                                items.append(listbox.get(idx))
                        plot_dtype_selections.append(items)
                    else:
                        plot_dtype_selections.append([listbox.get(idx) for idx in selected_indices])
                plot_dtype_prev_selected.clear()
                plot_dtype_prev_selected.extend(plot_dtype_selections)

            def numerator_listbox_on_select(event):
                plot_numerator_selections.clear()
                for i, listbox in enumerate(plot_numerator_listboxes):
                    selected_indices = listbox.curselection()
                    if len(plot_numerator_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in plot_numerator_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: plot_numerator_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in plot_numerator_prev_selected[i]:
                                items.append(listbox.get(idx))
                        plot_numerator_selections.append(items)
                    else:
                        plot_numerator_selections.append([listbox.get(idx) for idx in selected_indices])
                plot_numerator_prev_selected.clear()
                plot_numerator_prev_selected.extend(plot_numerator_selections)

            def denominator_listbox_on_select(event):
                plot_denominator_selections.clear()
                for i, listbox in enumerate(plot_denominator_listboxes):
                    selected_indices = listbox.curselection()
                    if len(plot_denominator_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in plot_denominator_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: plot_denominator_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in plot_denominator_prev_selected[i]:
                                items.append(listbox.get(idx))
                        plot_denominator_selections.append(items)
                    else:
                        plot_denominator_selections.append([listbox.get(idx) for idx in selected_indices])
                plot_denominator_prev_selected.clear()
                plot_denominator_prev_selected.extend(plot_denominator_selections)

            def condition_listbox_on_select(event):
                plot_condition_selections.clear()
                for i, listbox in enumerate(plot_condition_listboxes):
                    selected_indices = listbox.curselection()
                    if len(plot_condition_prev_selected) != 0:
                        items = []
                        for idx in selected_indices:
                            if listbox.get(idx) in plot_condition_prev_selected[i]:
                                items.append(listbox.get(idx))
                        items = sorted(items, key=lambda x: plot_condition_prev_selected[i].index(x))
                        for idx in selected_indices:
                            if listbox.get(idx) not in plot_condition_prev_selected[i]:
                                items.append(listbox.get(idx))
                        plot_condition_selections.append(items)
                    else:
                        plot_condition_selections.append([listbox.get(idx) for idx in selected_indices])
                plot_condition_prev_selected.clear()
                plot_condition_prev_selected.extend(plot_condition_selections)

            #def destroy_windows():
            #    for file_name in os.listdir(os.getcwd()):
            #        if file_name[0:4] == "temp":
            #            os.remove(file_name)
            #    root2.destroy()

            root2 = tk.Tk.__init__(self, *args, **kwargs)
            self.frame = VerticalScrolledFrame(root2)
            self.frame.pack()

            #background_image = ttk.PhotoImage(self.frame.interior, file="./GUI/resized_image.png")
            #bg = ttk.Label(self.frame.interior, image = background_image)
            #bg.place(x = 0,y = 0)

            # PLOT TYPE OPTIONS
            plot_type_selections = []
            plot_type_labels = []
            plot_type_option_menus = []
            plot_type_options = ["line", "two-axis", "jitter", "grouped jitter", "heatmap"]
            plot_type_frm = ttk.Frame(self.frame.interior)
            plot_type_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                plot_type_selection = ttk.StringVar()
                plot_type_selection.set("")
                plot_type_selections.append(plot_type_selection)
                plot_type_label = ttk.Label(plot_type_frm, text=f"Plot {i} type:")
                plot_type_label.pack(side=LEFT, anchor="w", padx=5)
                plot_type_option_menu = ttk.OptionMenu(plot_type_frm,
                                                       plot_type_selections[i-1],
                                                       *plot_type_options,
                                                       command=lambda value, index = i: plot_type_select(value, index))
                plot_type_option_menu.pack(side=LEFT, after=plot_type_label, anchor="w", padx=5)
                plot_type_option_menus.append(plot_type_option_menu)

            # PLOT DTYPE OPTIONS

            plot_dtype_listboxes = []
            plot_dtype_prev_selected =  []
            plot_dtype_selections = []
            plot_dtype_options = ["lum", "OD", "RLU", "BF_imaging", "CFP_imaging", "YFP_imaging", "texas_red_imaging",
                                  "CY5_imaging", "YFP", "CY5"]
            plot_dtype_options_lens = [len(x) for x in plot_dtype_options]
            max_len_dtype_ind = np.argmax(plot_dtype_options_lens)
            plot_dtype_frm = ttk.Frame(self.frame.interior)
            plot_dtype_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                # scrollbar = ttk.Scrollbar(frm2, orient=ttk.VERTICAL)
                plot_dtype_label = ttk.Label(plot_dtype_frm, text=f"Plot {i} data type:")
                plot_dtype_label.pack(side=LEFT, anchor="w")
                plot_dtype_listbox = ttk.Listbox(plot_dtype_frm,
                                                 selectmode=ttk.MULTIPLE,
                                                 height=len(plot_dtype_options),
                                                 width = plot_dtype_options_lens[max_len_dtype_ind],
                                                 exportselection=False)
                for option in plot_dtype_options:
                    plot_dtype_listbox.insert(ttk.END, option)
                plot_dtype_listbox.bind("<<ListboxSelect>>", dtype_listbox_on_select)
                plot_dtype_listbox.pack(side=ttk.LEFT, after = plot_dtype_label, padx=10)
                plot_dtype_listboxes.append(plot_dtype_listbox)

            # PLOT NUMERATOR OPTIONS

            plot_numerator_listboxes = []
            plot_numerator_prev_selected =  []
            plot_numerator_selections = []
            plot_numerator_options = ["lum", "OD", "RLU", "BF_imaging", "CFP_imaging", "YFP_imaging", "texas_red_imaging",
                                  "CY5_imaging", "YFP", "CY5"]
            plot_numerator_options_lens = [len(x) for x in plot_numerator_options]
            max_len_numerator_ind = np.argmax(plot_numerator_options_lens)
            plot_numerator_frm = ttk.Frame(self.frame.interior)
            plot_numerator_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                # scrollbar = ttk.Scrollbar(frm2, orient=ttk.VERTICAL)
                plot_numerator_label = ttk.Label(plot_numerator_frm, text=f"Plot {i} numerator(s):")
                plot_numerator_label.pack(side=LEFT, anchor="w")
                plot_numerator_listbox = ttk.Listbox(plot_numerator_frm,
                                                 selectmode=ttk.MULTIPLE,
                                                 height=len(plot_numerator_options),
                                                 width = plot_numerator_options_lens[max_len_numerator_ind],
                                                 exportselection=False)
                for option in plot_numerator_options:
                    plot_numerator_listbox.insert(ttk.END, option)
                plot_numerator_listbox.bind("<<ListboxSelect>>", numerator_listbox_on_select)
                plot_numerator_listbox.pack(side=ttk.LEFT, after = plot_numerator_label, padx=10)
                plot_numerator_listboxes.append(plot_numerator_listbox)

            # PLOT DENOMINATOR OPTIONS

            plot_denominator_listboxes = []
            plot_denominator_prev_selected =  []
            plot_denominator_selections = []
            plot_denominator_options = ["lum", "OD", "RLU", "BF_imaging", "CFP_imaging", "YFP_imaging", "texas_red_imaging",
                                  "CY5_imaging", "YFP", "CY5"]
            plot_denominator_options_lens = [len(x) for x in plot_denominator_options]
            max_len_denominator_ind = np.argmax(plot_denominator_options_lens)
            plot_denominator_frm = ttk.Frame(self.frame.interior)
            plot_denominator_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                # scrollbar = ttk.Scrollbar(frm2, orient=ttk.VERTICAL)
                plot_denominator_label = ttk.Label(plot_denominator_frm, text=f"Plot {i} denominator(s):")
                plot_denominator_label.pack(side=LEFT, anchor="w")
                plot_denominator_listbox = ttk.Listbox(plot_denominator_frm,
                                                 selectmode=ttk.MULTIPLE,
                                                 height=len(plot_denominator_options),
                                                 width = plot_denominator_options_lens[max_len_denominator_ind],
                                                 exportselection=False)
                for option in plot_denominator_options:
                    plot_denominator_listbox.insert(ttk.END, option)
                plot_denominator_listbox.bind("<<ListboxSelect>>", denominator_listbox_on_select)
                plot_denominator_listbox.pack(side=ttk.LEFT, after = plot_denominator_label, padx=10)
                plot_denominator_listboxes.append(plot_denominator_listbox)


            # PLOT CONDITION OPTIONS
            plot_condition_listboxes = []
            plot_condition_selections = []
            plot_condition_prev_selected = []
            plot_condition_options_lens = [len(x) for x in condition_names]
            max_len_condition_ind = np.argmax(plot_condition_options_lens)
            plot_condition_frm = ttk.Frame(self.frame.interior)
            plot_condition_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                plot_condition_label = ttk.Label(plot_condition_frm, text=f"Plot {i} condition(s):")
                plot_condition_label.pack(side=LEFT, anchor="w", padx=5)
                plot_condition_listbox = ttk.Listbox(plot_condition_frm,
                                                     selectmode=ttk.MULTIPLE,
                                                     height=len(condition_names),
                                                     width = plot_condition_options_lens[max_len_condition_ind]*2,
                                                     exportselection=False)
                for option in condition_names:
                    plot_condition_listbox.insert(ttk.END, option)
                plot_condition_listbox.bind("<<ListboxSelect>>", condition_listbox_on_select)

                plot_condition_listboxes.append(plot_condition_listbox)
                plot_condition_listbox.pack(side=ttk.LEFT, after = plot_condition_label, padx=10)

            # PLOT NORMALIZATION OPTIONS
            plot_normalization_selections = []
            plot_normalization_labels = []
            plot_normalization_option_menus = []
            plot_normalization_frm = ttk.Frame(self.frame.interior)
            plot_normalization_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                plot_normalization_selection = ttk.StringVar()
                plot_normalization_selection.set("")
                plot_normalization_selections.append(plot_normalization_selection)
                plot_normalization_label = ttk.Label(plot_normalization_frm, text=f"Plot {i} normalization:")
                plot_normalization_label.pack(side=LEFT, anchor="w", padx=5)
                plot_normalization_option_menu = ttk.OptionMenu(plot_normalization_frm,
                                                                plot_normalization_selections[i-1],
                                                                *condition_names+[""],
                                                                command=lambda value, index=i: plot_normalization_select(value,index))
                plot_normalization_option_menu.pack(side=LEFT, after=plot_normalization_label, anchor="w", padx=10)
                plot_normalization_option_menus.append(plot_normalization_option_menu)

            # NORMALIZATION METHOD OPTIONS
            normalization_method_listboxes = []
            normalization_method_selections = []
            normalization_method_prev_selected =  []
            normalization_method_options = ["", "fold-change", "percent"]
            normalization_method_options_lens = [len(x) for x in normalization_method_options]
            max_len_norm_ind = np.argmax(normalization_method_options_lens)
            normalization_method_frm = ttk.Frame(self.frame.interior)
            normalization_method_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                # scrollbar = ttk.Scrollbar(frm2, orient=ttk.VERTICAL)
                normalization_method_label = ttk.Label(normalization_method_frm, text=f"Plot {i} normalization method(s):")
                normalization_method_label.pack(side=LEFT, anchor="w")
                normalization_method_listbox = ttk.Listbox(normalization_method_frm,
                                                 selectmode=ttk.MULTIPLE,
                                                 height=len(normalization_method_options),
                                                 width = normalization_method_options_lens[max_len_norm_ind],
                                                 exportselection=False)
                for option in normalization_method_options:
                    normalization_method_listbox.insert(ttk.END, option)
                normalization_method_listbox.bind("<<ListboxSelect>>", normalization_method_listbox_on_select)
                normalization_method_listbox.pack(side=ttk.LEFT, after = normalization_method_label, padx=10)
                normalization_method_listboxes.append(normalization_method_listbox)

            # PLOT XAXIS OPTIONS
            plot_xaxis_selections = []
            plot_xaxis_labels = []
            plot_xaxis_option_menus = []
            plot_xaxis_options = ["default", "Time", "OD"]
            plot_xaxis_frm = ttk.Frame(self.frame.interior)
            plot_xaxis_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                plot_xaxis_selection = ttk.StringVar()
                plot_xaxis_selection.set("")
                plot_xaxis_selections.append(plot_xaxis_selection)
                plot_xaxis_label = ttk.Label(plot_xaxis_frm, text=f"Plot {i} x-axis:")
                plot_xaxis_label.pack(side=LEFT, anchor="w", padx=5)
                plot_xaxis_option_menu = ttk.OptionMenu(plot_xaxis_frm,
                                                       plot_xaxis_selections[i-1],
                                                       *plot_xaxis_options,
                                                       command=lambda value, index = i: plot_xaxis_select(value, index))
                plot_xaxis_option_menu.pack(side=LEFT, after=plot_xaxis_label, anchor="w", padx=5)
                plot_xaxis_option_menus.append(plot_xaxis_option_menu)

            # PLOT SCALE OPTIONS
            plot_scale_listboxes = []
            plot_scale_selections = []
            plot_scale_prev_selected =  []
            plot_scale_options = ["log", "linear"]
            plot_scale_options_lens = [len(x) for x in plot_scale_options]
            plot_scale_ind = np.argmax(plot_scale_options_lens)
            plot_scale_frm = ttk.Frame(self.frame.interior)
            plot_scale_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                # scrollbar = ttk.Scrollbar(frm2, orient=ttk.VERTICAL)
                plot_scale_label = ttk.Label(plot_scale_frm, text=f"Plot {i} y-scale:")
                plot_scale_label.pack(side=LEFT, anchor="w")
                plot_scale_listbox = ttk.Listbox(plot_scale_frm,
                                                 selectmode=ttk.MULTIPLE,
                                                 height=len(plot_scale_options),
                                                 width = plot_scale_options_lens[plot_scale_ind],
                                                 exportselection=False)
                for option in plot_scale_options:
                    plot_scale_listbox.insert(ttk.END, option)
                plot_scale_listbox.bind("<<ListboxSelect>>", plot_scale_listbox_on_select)
                plot_scale_listbox.pack(side=ttk.LEFT, after = plot_scale_label, padx=10)
                plot_scale_listboxes.append(plot_scale_listbox)

            # DEFAULT PLOT COLOR
            plot_colors_frm = ttk.Frame(self.frame.interior)
            plot_colors_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            plot_color_entries = []
            for i in range(1, num_plots + 1):
                plot_colors_label = ttk.Label(plot_colors_frm,
                                             text=f"Plot {i} default color (hex):")
                plot_colors_label.pack(side=LEFT, anchor="w", padx=5)
                plot_colors_entry = ttk.Entry(plot_colors_frm)
                plot_colors_entry.pack(side=LEFT, after=plot_colors_label, anchor="w", padx=5)
                plot_color_entries.append(plot_colors_entry)

            # PLOT FONT OPTIONS
            plot_font_selections = []
            plot_font_labels = []
            plot_font_option_menus = []
            plot_font_options = ["Helvetica", "Computer modern"]
            plot_font_frm = ttk.Frame(self.frame.interior)
            plot_font_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            for i in range(1, num_plots + 1):
                plot_font_selection = ttk.StringVar()
                plot_font_selection.set("")
                plot_font_selections.append(plot_font_selection)
                plot_font_label = ttk.Label(plot_font_frm, text=f"Plot {i} font:")
                plot_font_label.pack(side=LEFT, anchor="w", padx=5)
                plot_font_option_menu = ttk.OptionMenu(plot_font_frm,
                                                       plot_font_selections[i-1],
                                                       *plot_font_options,
                                                       command=lambda value, index = i: plot_font_select(value, index))
                plot_font_option_menu.pack(side=LEFT, after=plot_font_label, anchor="w", padx=5)
                plot_font_option_menus.append(plot_font_option_menu)

            # PLOT TITLES
            plot_title_frm = ttk.Frame(self.frame.interior)
            plot_title_frm.pack(side=TOP, anchor="w", padx=5, pady=5)

            plot_title_entries = []
            for i in range(1, num_plots + 1):
                plot_title_label = ttk.Label(plot_title_frm,
                                             text=f"Plot {i} title:")
                plot_title_label.pack(side=LEFT, anchor="w", padx=5)
                plot_title_entry = ttk.Entry(plot_title_frm)
                plot_title_entry.pack(side=LEFT, after=plot_title_label, anchor="w", padx=5)
                plot_title_entries.append(plot_title_entry)

            # PLOT YLABS
            plot_ylabs_frm = ttk.Frame(self.frame.interior)
            plot_ylabs_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)

            plot_ylab_entries = []
            for i in range(1, num_plots + 1):
                plot_ylabs_label = ttk.Label(plot_ylabs_frm,
                                             text=f"Plot {i} y-label:")
                plot_ylabs_label.pack(side=LEFT, anchor="w", padx=5)
                plot_ylabs_entry = ttk.Entry(plot_ylabs_frm)
                plot_ylabs_entry.pack(side=LEFT, after=plot_ylabs_label, anchor="w", padx=5)
                plot_ylab_entries.append(plot_ylabs_entry)

            # PLOT XLABS
            plot_xlabs_frm = ttk.Frame(self.frame.interior)
            plot_xlabs_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            plot_xlab_entries = []
            for i in range(1, num_plots + 1):
                plot_xlabs_label = ttk.Label(plot_xlabs_frm,
                                             text=f"Plot {i} x-label:")
                plot_xlabs_label.pack(side=LEFT, anchor="w", padx=5)
                plot_xlabs_entry = ttk.Entry(plot_xlabs_frm)
                plot_xlabs_entry.pack(side=LEFT, after=plot_xlabs_label, anchor="w", padx=5)
                plot_xlab_entries.append(plot_xlabs_entry)

            # PLOT YTICKS
            plot_yticks_frm = ttk.Frame(self.frame.interior)
            plot_yticks_frm.pack(side=TOP, anchor = "w", padx=5, pady=5)

            plot_ytick_entries = []
            for i in range(1, num_plots + 1):
                plot_yticks_label = ttk.Label(plot_yticks_frm,
                                             text=f"Plot {i} y-ticks:")
                plot_yticks_label.pack(side=LEFT, anchor="w", padx=5)
                plot_yticks_entry = ttk.Entry(plot_yticks_frm)
                plot_yticks_entry.pack(side=LEFT, after=plot_yticks_label, anchor="w", padx=5)
                plot_ytick_entries.append(plot_yticks_entry)

            # PLOT XTICKS
            plot_xticks_frm = ttk.Frame(self.frame.interior)
            plot_xticks_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            plot_xtick_entries = []
            for i in range(1, num_plots + 1):
                plot_xticks_label = ttk.Label(plot_xticks_frm,
                                             text=f"Plot {i} x-ticks:")
                plot_xticks_label.pack(side=LEFT, anchor="w", padx=5)
                plot_xticks_entry = ttk.Entry(plot_xticks_frm)
                plot_xticks_entry.pack(side=LEFT, after=plot_xticks_label, anchor="w", padx=5)
                plot_xtick_entries.append(plot_xticks_entry)

            # COLOR LABELS
            color_label_frm = ttk.Frame(self.frame.interior)
            color_label_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            color_label_entries = []
            for i in range(1, num_plots + 1):
                color_label_label = ttk.Label(color_label_frm,
                                                text=f"Plot {i} colorbar label:")
                color_label_label.pack(side=LEFT, anchor="w", padx=5)
                color_label_entry = ttk.Entry(color_label_frm)
                color_label_entry.pack(side=LEFT, after=color_label_label, anchor="w", padx=5)
                color_label_entries.append(color_label_entry)

            # PLOT SIZES
            plot_sizes_frm = ttk.Frame(self.frame.interior)
            plot_sizes_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            plot_size_entries = []
            for i in range(1, num_plots + 1):
                plot_sizes_label = ttk.Label(plot_sizes_frm,
                                             text=f"Plot {i} size:")
                plot_sizes_label.pack(side=LEFT, anchor="w", padx=5)
                plot_sizes_entry = ttk.Entry(plot_sizes_frm)
                plot_sizes_entry.pack(side=LEFT, after=plot_sizes_label, anchor="w", padx=5)
                plot_size_entries.append(plot_sizes_entry)

            # PLOT YLIMS
            ylims_frm = ttk.Frame(self.frame.interior)
            ylims_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            ylim_entries = []
            for i in range(1, num_plots + 1):
                ylims_label = ttk.Label(ylims_frm,
                                             text=f"Plot {i} ylimits (jitter only):")
                ylims_label.pack(side=LEFT, anchor="w", padx=5)
                ylims_entry = ttk.Entry(ylims_frm)
                ylims_entry.pack(side=LEFT, after=ylims_label, anchor="w", padx=5)
                ylim_entries.append(ylims_entry)

            # DOSE CONCS
            dose_concs_frm = ttk.Frame(self.frame.interior)
            dose_concs_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            dose_conc_entries = []
            for i in range(1, num_plots + 1):
                dose_concs_label = ttk.Label(dose_concs_frm,
                                             text=f"Plot {i} dose concentrations:")
                dose_concs_label.pack(side=LEFT, anchor="w", padx=5)
                dose_concs_entry = ttk.Entry(dose_concs_frm)
                dose_concs_entry.pack(side=LEFT, after=dose_concs_label, anchor="w", padx=5)
                dose_conc_entries.append(dose_concs_entry)

            # PLOT FILENAMES
            plot_filename_frm = ttk.Frame(self.frame.interior)
            plot_filename_frm.pack(side=TOP, anchor="w", padx=5, pady=5)
            plot_filename_entries = []
            for i in range(1, num_plots + 1):
                plot_filename_label = ttk.Label(plot_filename_frm,
                                                text=f"Plot {i} filename:")
                plot_filename_label.pack(side=LEFT, anchor="w", padx=5)
                plot_filename_entry = ttk.Entry(plot_filename_frm)
                plot_filename_entry.pack(side=LEFT, after=plot_filename_label, anchor="w", padx=5)
                plot_filename_entries.append(plot_filename_entry)

            # SAVE TO JSON FILE
            save_to_json_btn = ttk.Button(self.frame.interior,
                                          text="Save to JSON",
                                          command=save_to_json)
            save_to_json_btn.pack(side=TOP, anchor = "s", padx = 10)
            # EXIT
            #quit_btn = ttk.Button(self.frame.interior, text="Done",
            #                      command=destroy_windows)
            #quit_btn.pack(side=TOP, padx = 10)

    app = PlotOptions()
    app.mainloop()
