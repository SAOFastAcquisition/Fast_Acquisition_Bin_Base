import tkinter as tk
from PIL import Image as PilImage
from PIL import ImageTk


class Wind:

    def __init__(self, width=300, height=200, title='Main Window', resize=(True, True), icon=None):
        self.wind = tk.Tk()
        self.wind.geometry(f'{width}x{height}+200+200')
        self.wind.title(title)
        self.wind.resizable(resize[0], resize[1])

        if icon:
            photo_icon = tk.PhotoImage(file=icon)
            self.wind.iconphoto(False, photo_icon)

        self.lbl = tk.Label(self.wind, text='Attention!', bg='#245e11', fg='#2a2d73', relief='groove',
                            wraplength=55, font='Consoles 16')
        self.btn1 = tk.Button(self.wind, text='Initial Button', width=15, height=1, bg='blue', fg='yellow',
                              relief='groove', activebackground='red', font=('Consolas', 12, 'bold'))
        text_var = tk.StringVar(value='Name Var Button')
        self.btn2 = tk.Button(self.wind, textvariable=text_var, width=15, height=1, bg='blue', fg='yellow',
                              relief='groove', activebackground='red', font=('Consolas', 12, 'bold'))
        self.image1 = tk.PhotoImage(file='ghost.gif')
        self.lbl_pic = tk.Label(self.wind, text='ghost', image=self.image1, compound='left')
        # *********************** Кнопка с картинкой *************************
        img1 = PilImage.open('ghost.gif')                # Считываем изображение
        img1 = img1.resize((20, 20), PilImage.ANTIALIAS) # Приводим его в размер, пригодный для кнопки
        self.photo_img1 = ImageTk.PhotoImage(img1)       # Формируем переменную для встраивания в кнопку
        self.btn3 = tk.Button(self.wind, text='ghost', width=50, height=25, image=self.photo_img1, compound='left',
                              command=self.btn3_action)
        # *****************************************************************
        # ********************** Кнопки выхода из окна ********************
        self.destr_btn = tk.Button(self.wind, text='Destroy Window', bg='red', command=self.wind.destroy)
        self.quit_btn = tk.Button(self.wind, text='Quit', bg='purple', command=quit)
        # *****************************************************************
    def run(self):
        self.show_widget()
        # self.draw_widgets()
        self.wind.mainloop()

    def btn3_action(self):
        tk.Label(self.wind, width=15, height=2, bg='red', text='First').pack(side='left', padx=10)
        self.lbl.config(text='Button pushed!', bg='purple')
        print(self.btn3)
        pass

    def create_child(self, width, height, title='Child Window', resize=(False, False), icon=None):
        Child_Wind(self.wind, width, height, title, resize, icon)

    def show_widget(self):
        self.lbl.pack(anchor='se', padx=100, pady=10)
        self.btn1.pack()
        self.btn2.pack()
        self.btn3.pack()
        self.destr_btn.pack()
        self.quit_btn.pack()

    def draw_widgets_pack(self):
        frame1 = tk.LabelFrame(self.wind, text='Top Frame', padx=0, pady=0)
        frame2 = tk.Frame(self.wind, padx=15, pady=10)
        frame1.pack(side='top', fill='both')
        frame2.pack(ipadx=200, ipady=50, fill='both', anchor='se')

        tk.Label(frame1, width=15, height=2, bg='red', text='First').pack(side='left', padx=10)
        tk.Label(frame1, width=15, height=2, bg='yellow', text='Second').pack(side='right')
        tk.Label(frame2, width=15, height=2, bg='green', text='Third').pack(side='left', pady=20)
        tk.Label(frame2, width=15, height=2, bg='cyan', text='Fourth').pack(side='left')

    def draw_widgets_place(self):
        frame1 = tk.Frame(self.wind, borderwidth=20, relief='sunken', width=300, height=100)
        frame2 = tk.Frame(self.wind, borderwidth=20, relief='sunken', width=300, height=100)
        frame1.place(x=0, y=0)
        frame2.place(x=0, y=200)

        # tk.Label(self.wind, width=15, height=2, bg='red', text='First').place(x=10, y=10, relwidth=0.35, relheight=0.1)
        # tk.Label(self.wind, width=15, height=2, bg='yellow', text='Second').place(relx=0.5, rely=0.1, width=75, height=30)
        tk.Label(frame1, width=15, height=2, bg='green', text='Third').place(x=10, y=10, relwidth=0.35, relheight=0.35)
        tk.Label(frame2, width=15, height=2, bg='cyan', text='Fourth').place(x=10, y=10, relwidth=0.35, relheight=0.35)

    def draw_widgets(self):
        # frame1 = tk.Frame(self.wind, borderwidth=20, relief='sunken', width=300, height=100)
        # frame2 = tk.Frame(self.wind, borderwidth=20, relief='sunken', width=300, height=100)
        # frame1.place(x=0, y=0)
        # frame2.place(x=0, y=200)

        tk.Label(self.wind, width=15, height=2, bg='red', text='First').grid(row=0, column=0)
        tk.Label(self.wind, width=15, height=2, bg='yellow', text='Second').grid(row=1, column=1)
        # tk.Label(frame1, width=15, height=2, bg='green', text='Third').place(x=10, y=10, relwidth=0.35, relheight=0.35)
        # tk.Label(frame2, width=15, height=2, bg='cyan', text='Fourth').place(x=10, y=10, relwidth=0.35, relheight=0.35)


class Child_Wind:

    def __init__(self, parent, width=300, height=200, title='Main Window', resize=(False, False), icon=None):
        self.ch_wind = tk.Toplevel(parent)
        self.ch_wind.geometry(f'{width}x{height}+500+100')
        self.ch_wind.title(title)
        self.ch_wind.resizable(resize[0], resize[1])
        if icon:
            self.ch_wind.iconbitmap(icon)

        self.grab_focus()

    def grab_focus(self):
        self.ch_wind.grab_set()
        self.ch_wind.focus_set()
        self.ch_wind.wait_window()


if __name__ == '__main__':
    win = Wind(500, 400, icon='ghost.gif')
    # win.create_child(300, 200, icon='icon.ico')
    win.run()
    pass
