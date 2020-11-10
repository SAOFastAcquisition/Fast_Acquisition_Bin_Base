import tkinter as tk
from tkinter import messagebox as mb
from PIL import Image as pilImage
from PIL import ImageTk, ImageOps


class Wind:

    def __init__(self, width, height, title='Main Window', resize=(False, False), icon=None, ):
        self.wind = tk.Tk()
        self.wind.geometry(f'{width}x{height}+200+200')
        self.wind.title(title)
        self.wind.resizable(resize[0], resize[1])
        if icon:
            photo_icon = tk.PhotoImage(file=icon)
            self.wind.iconphoto(False, photo_icon)

    #     Radiobuttons
        img = pilImage.open('ghost.gif')
        neg = ImageOps.invert(img.convert('RGB'))
        img = img.resize((20,20), pilImage.ANTIALIAS)
        neg = neg.resize((20,20), pilImage.ANTIALIAS)
        self.ghost = ImageTk.PhotoImage(img)
        self.neg_ghost = ImageTk.PhotoImage(neg)

        self.choice = tk.IntVar(value=1)
        self.rbtn1 = tk.Radiobutton(self.wind, text='One', variable=self.choice, value=1, width=70, height=10,
                                    bg='green', relief='groove', bd=5, image=self.ghost, selectimage=self.neg_ghost,
                                    justify='right')
        self.rbtn2 = tk.Radiobutton(self.wind, text='Two', variable=self.choice, value=2, width=10, height=3,
                                    bg='blue', relief='groove', bd=5)
        self.rbtn3 = tk.Radiobutton(self.wind, text='Three', variable=self.choice, value=3, width=10, height=3,
                                    bg='cyan', relief='groove', bd=5)

        self.choice1 = tk.IntVar(value=31)
        self.rbtn11 = tk.Radiobutton(self.wind, text='One', variable=self.choice1, value=11, width=10, height=3,
                                    bg='green', relief='sunken', bd=5)
        self.rbtn21 = tk.Radiobutton(self.wind, text='Two', variable=self.choice1, value=21, width=10, height=3,
                                    bg='green', relief='sunken', bd=5)
        self.rbtn31 = tk.Radiobutton(self.wind, text='Three', variable=self.choice1, value=31, width=10, height=3,
                                    bg='green', relief='sunken', bd=5)
        self.rbtn41 = tk.Radiobutton(self.wind, text='Four', variable=self.choice1, value=41, width=10, height=3,
                                    command=self.command41, bg='green', relief='sunken', bd=5)

        self.choice2 = tk.IntVar(value=12)
        self.colour1 = 'red'
        self.colour2 = 'purple'
        self.colour3 = 'grey'
        self.rbtn12 = tk.Radiobutton(self.wind, text='One', variable=self.choice2, value=12, width=70, height=10,
                                    bg='green', relief='groove', bd=5, image=self.ghost, selectimage=self.neg_ghost,
                                    command=self.change_colour, justify='right')
        self.rbtn22 = tk.Radiobutton(self.wind, text='Two', variable=self.choice2, value=22, width=10, height=3,
                                    bg='blue', relief='groove', bd=5, command=self.change_colour)
        self.rbtn32 = tk.Radiobutton(self.wind, text='Three', variable=self.choice2, value=32, width=10, height=3,
                                    bg='cyan', relief='groove', bd=5, command=self.change_colour)

        self.rbtn4 = tk.Radiobutton(self.wind)

    # Buttons
        self.btn_del = tk.Button(self.wind, text='Quit', command=quit)
        self.btn_save = tk.Button(self.wind, text='Save rbtn41', command=self.command41)

    def widget_draw(self):
        self.rbtn1.grid(column=0, row=0)
        self.rbtn2.grid(column=0, row=1)
        self.rbtn3.grid(column=0, row=2)

        self.rbtn11.grid(column=1, row=0)
        self.rbtn21.grid(column=1, row=1)
        self.rbtn31.grid(column=1, row=2)
        self.rbtn41.grid(column=1, row=3)

        self.rbtn12.grid(column=2, row=0)
        self.rbtn22.grid(column=2, row=1)
        self.rbtn32.grid(column=2, row=2)

        self.btn_del.grid(column=0, row=5, columnspan=2)
        self.btn_save.grid(column=0, row=3)
        self.wind.mainloop()

    def del_entry(self):
        self.entr1.delete(0, tk.END)

    def command41(self):
        choice= self.choice1.get()
        mb.showinfo('Info', f'Choice = {choice}')

    def change_colour(self):
        choice2 = self.choice2.get()
        if choice2 == 12:
            self.rbtn12.config(bg=self.colour1)
        else:
            self.rbtn12.config(bg='green')
        if choice2 == 22:
            self.rbtn22.config(bg=self.colour2)
        else:
            self.rbtn22.config(bg='blue')
        if choice2 == 32:
            self.rbtn32.config(bg=self.colour3)
        else:
            self.rbtn32.config(bg='cyan')


if __name__ == '__main__':
    wndw = Wind(500, 400, icon='ghost.gif')
    wndw.widget_draw()