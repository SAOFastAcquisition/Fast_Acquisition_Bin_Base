import tkinter as tk
from tkinter import messagebox as mb
from PIL import Image as pilImage
from PIL import ImageTk, ImageOps


class Wind:

    def __init__(self, width, height, title='Calculator', resize=(False, False), icon=None, ):
        self.wind = tk.Tk()
        self.wind.geometry(f'{width}x{height}+200+400')
        self.wind['bg'] = 'cyan'
        self.wind.title(title)
        self.wind.resizable(resize[0], resize[1])
        if icon:
            photo_icon = tk.PhotoImage(file=icon)
            self.wind.iconphoto(False, photo_icon)

        self.calc = tk.Entry(self.wind, width=25, borderwidth=5, justify='right')#.grid(column=0, row=0, columnspan=4)

        # Buttons
        for i in range(10):
            if i == 0:
                tk.Button(self.wind, text=f'{i}', font=('consolas', 20), width=3, height=1, relief='groove',
                                      command=lambda: self.ins_num(i)).grid(row=7 - (i + 2) // 3, column=i, padx=5, pady=5)
            else:
                tk.Button(self.wind, text=f'{i}', font=('consolas', 20), width=3, height=1, relief='groove'). \
                    grid(row=7 - (i + 2) // 3, column=(i - 1) % 3, padx=5, pady=5)

        self.btn_del = tk.Button(self.wind, text='Quit', command=quit)

    def ins_num(self, i):
        num = self.calc.get()
        num = num + str(i)
        self.calc.delete(0, tk.END)
        self.calc.insert(0, num)

    def widget_draw(self):
        self.calc.grid(column=0, row=0, columnspan=4)
        # Buttons

        # self.btn_del.grid(column=0, row=5, columnspan=2)
        self.wind.mainloop()

    def del_entry(self):
        self.entr1.delete(0, tk.END)


if __name__ == '__main__':
    wndw = Wind(500, 400, icon='ghost.gif')
    wndw.widget_draw()