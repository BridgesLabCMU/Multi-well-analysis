import tkinter as ttk
from tkinter import *

sample_cells = set()

def num_to_letter(num):
    if num == 0:
        return "A"
    elif num == 1:
        return "B"
    elif num == 2:
        return "C"
    elif num == 3:
        return "D"
    elif num == 4:
        return "E"
    elif num == 5:
        return "F"
    elif num == 6:
        return "G"
    elif num == 7:
        return "H"

class Table(ttk.Frame):
    def __init__(self, master, rows, columns):
        super().__init__(master)
        self.rows = rows
        self.columns = columns
        self.cell_width = 50
        self.cell_height = 30
        self.selected_cells = set()
        self.create_widgets()

    def create_widgets(self):
        self.canvas = ttk.Canvas(self, borderwidth=0, highlightthickness=0)
        self.canvas.pack(side="top", fill="both", expand=True)

        self.canvas.bind("<Button-1>", self.on_click)

        self.draw_table()

    def draw_table(self):
        for row in range(self.rows):
            for col in range(self.columns):
                x0 = col * self.cell_width
                y0 = row * self.cell_height
                x1 = x0 + self.cell_width
                y1 = y0 + self.cell_height

                cell_id = self.canvas.create_rectangle(x0, y0, x1, y1, outline="black", fill="white")
                self.canvas.tag_bind(cell_id, "<Button-1>", lambda event, cell=(row, col): self.cell_clicked(cell))

    def cell_clicked(self, cell):
        row, col = cell
        
        r = num_to_letter(row)
        c = str(col + 1)
        cell_val = r + c
        
        cell_id = self.get_cell_id(row, col)
        if cell_id not in self.selected_cells:
            self.selected_cells.add(cell_id)
            self.canvas.itemconfig(cell_id, fill="lightblue")
            sample_cells.add(cell_val)
        else:
            self.selected_cells.remove(cell_id)
            self.canvas.itemconfig(cell_id, fill="white")
            sample_cells.remove(cell_val)

        
    def on_click(self, event):
        self.start_cell = self.get_coords(event)
        self.end_cell = self.get_coords(event) + (self.cell_width,self.cell_height)
        
    
    def get_coords(self, event):
        col = event.x // self.cell_width
        row = event.y // self.cell_height
        return row, col

    def get_cell_id(self, row, col):
        return (row * self.columns) + col + 1

    def clear_color(self):
        for cell in self.selected_cells:
            self.canvas.itemconfig(cell, fill="white")
        self.selected_cells = set()

