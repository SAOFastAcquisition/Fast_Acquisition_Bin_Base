import tkinter as tk
from tkinter import messagebox as mb

class Wind():

    def __init__(self, width, height, title='Main Window', resize=(False, False), icon=None, ):
        self.wind = tk.Tk()
        self.wind.geometry(f'{width}x{height}+200+200')
        self.wind.title(title)
        self.wind.resizable(resize[0], resize[1])
        if icon:
            photo_icon = tk.PhotoImage(file=icon)
            self.wind.iconphoto(False, photo_icon)

        self.lbl1 = tk.Label(self.wind, text='Name', width=15, height=3, bg='yellow', font=('Consolas', 16, 'bold'))
        text_var = tk.StringVar(value='prefix')
        self.entr1 = tk.Entry(self.wind, justify='right', font=('curier', 12), bg='purple', fg='yellow',
                              textvariable=text_var)
        self.lbl2 = tk.Label(self.wind, text='Password', width=15, height=3, bg='blue', font=('Consolas', 16, 'bold'))
        self.entr2 = tk.Entry(self.wind, show='+')
        self.btn1 = tk.Button(self.wind, text='Get Entry', width=15, height=3, bg='green', command=self.get_entry)
        self.btn2 = tk.Button(self.wind, text='Delete Entry', width=15, height=3, bg='green', command=self.del_entry)


    def widget_draw(self):
        self.lbl1.grid(row=0, column=0, sticky='w')
        self.lbl2.grid(row=1, column=0, sticky='e')
        self.entr1.grid(row=0, column=1)
        self.entr2.grid(row=1, column=1)
        self.btn1.grid(row=2, column=0, columnspan=2, stick='we')
        self.btn2.grid(row=3, column=0, columnspan=2, stick='we')
        self.wind.mainloop()

    def get_entry(self):
        word = self.entr1.get()
        pw = self.entr2.get()
        if word:
            print(word)
            mb.showinfo('Entry:', word)
        else:
            print('Empty entry')
        if pw:
            print(pw)
        else:
            print('Empty entry')
        self.entr1.delete(0, tk.END)
        self.entr2.delete(0, tk.END)

    def del_entry(self):
        self.entr1.delete(0, tk.END)


if __name__ == '__main__':
    wndw = Wind(500, 400, icon='ghost.gif')
    wndw.widget_draw()