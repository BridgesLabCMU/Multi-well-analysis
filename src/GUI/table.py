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
def letter_to_num(letter):
    if letter == "A":
        return 0
    elif letter == "B":
        return 1
    elif letter == "C":
        return 2
    elif letter == "D":
        return 3
    elif letter == "E":
        return 4
    elif letter == "F":
        return 5
    elif letter == "G":
        return 6
    elif letter == "H":
        return 7

class Table(ttk.Frame):
    def __init__(self, master, rows, columns):
        super().__init__(master)
        self.rows = rows - 1
        self.columns = columns - 1
        self.cell_width = 50
        self.cell_height = 30
        self.selected_cells = set()
        self.create_widgets()

    def create_widgets(self):
        self.canvas = ttk.Canvas(self, borderwidth=0, highlightthickness=0)
        self.canvas.pack(side="top", fill="both", expand=True)
        self.draw_table()
        self.canvas.bind("<Button-1>", self.on_click)


    def draw_table(self):
        for row in range(self.rows):
            for col in range(self.columns):
                x0 = col * self.cell_width
                y0 = row * self.cell_height
                x1 = x0 + self.cell_width
                y1 = y0 + self.cell_height



                if col == 0:
                    if row >= 1:
                        letter = num_to_letter(row - 1)
                        # self.canvas.create_rectangle(x0, y0, x1, y1, outline="black", fill="lightgreen")
                        self.canvas.create_text(x0 + 25, y1 - 15, text=letter)
                    continue
                if row == 0:
                    if col >= 1:
                        # self.canvas.create_rectangle(x0, y0, x1, y1, outline="black", fill="lightgreen")
                        self.canvas.create_text(x1 - 25, y0 + 15, text = str(col))
                    continue
                self.canvas.create_rectangle(x0, y0, x1, y1, outline="black", fill="white")

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
        self.canvas.bind("<ButtonRelease-1>", lambda event, cell=(self.start_cell[0], self.end_cell[1]): self.on_release(cell, event))

    def on_release(self, cell, event):
        click_row, click_col = cell
        self.release_start_cell = self.get_coords(event)
        release_row, release_col = self.release_start_cell


        if click_row == 0 or click_col == 0 or release_row == 0 or release_col == 0:
            return
        if click_row > release_row:
            temp_row = click_row
            click_row = release_row
            release_row = temp_row
        if click_col > release_col:
            temp_col = click_col
            click_col = release_col
            release_col = temp_col


        for i in range(click_row, release_row + 1):
            for j in range(click_col - 1, release_col):
                cell_id = self.get_cell_id(i,j)
                r = num_to_letter(i - 1)
                c = str(j + 1)
                cell_val = (r  + c)
                if self.canvas.itemconfig(cell_id)["fill"][4] == "white":
                    if cell_id not in self.selected_cells:
                        self.selected_cells.add(cell_id)
                        self.canvas.itemconfig(cell_id, fill = "lightblue")
                        sample_cells.add(cell_val)
                elif self.canvas.itemconfig(cell_id)["fill"][4] == "lightblue":
                    self.selected_cells.remove(cell_id)
                    self.canvas.itemconfig(cell_id, fill="white")
                    sample_cells.remove(cell_val)

    def get_coords(self, event):
        col = event.x // self.cell_width
        row = event.y // self.cell_height
        return row, col

    def get_cell_id(self, row, col):
        return (row * self.columns) + col + 1

    def clear_color(self, cells):
        for cell in self.selected_cells:
            self.canvas.itemconfig(cell, fill="white")
        for id in range(13,110):
            if id%12 == 0:
                continue
            self.canvas.itemconfig(id, fill="white")
        self.selected_cells = set()

    def gray_out_cells(self, cells):
        for cell in cells:
            i = letter_to_num(cell[0]) + 1
            j = int(cell[1:]) - 1
            cell_id = self.get_cell_id(i, j)
            self.canvas.itemconfig(cell_id, fill = "darkgray")

