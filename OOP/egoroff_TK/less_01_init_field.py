
import tkinter as tk
from PIL import Image, ImageTk


def clicked():
    lbl1.config(text="I asked - ")


def clicked_entr1():
    res = "Input {}".format(input_txt.get())
    lbl2.config(text=res)


wind = tk.Tk()
wight = 600
height = 500
wind.title('Main Window')


photo = tk.PhotoImage(file="ghost.gif")
wind.iconphoto(False, photo)
# wind.iconbitmap('ghost.ico')
wind.config(bg='#00E600')
wind.geometry(f'{wight}x{height}+400+100')
wind.resizable(True, False)
wind.minsize(400, 500)
wind.maxsize(800, 700)

lbl1 = tk.Label(wind, text='Coeff', font=('TimesRoman', 16))
lbl2 = tk.Label(wind, text='Coeff1', font=('TimesRoman', 16))
lbl1.grid(column=0, row=1)
lbl2.grid(column=2, row=1)

input_txt = tk.Entry(wind, width=16)
input_txt.focus
input_txt.grid(column=0, row=0)
# lbl.pack()


def clicked_entr():

    res = "Input {}".format(input_txt.get())
    input_txt.config(state='disabled')
    lbl2.config(text=res)


def click_counter():
    global count
    count += 1
    btn4['text'] = f'Counter {count}'
    tk.Label(wind, text=f'{count}', font=('Times Roman', 16)).grid(row=4, column=1+count)



btn1 = tk.Button(wind, text="Don't push!!!", font=('Times Roman', 16),
                 command=clicked_entr
                 )
btn1.grid(column=1, row=0)

btn2 = tk.Button(wind, text="Don't push!!!", font=('Times Roman', 16),
                 command=clicked
                 )
btn2.grid(column=1, row=1)

btn3 = tk.Button(wind, text='Create new label', font=('Times Roman', 16),
                 command=lambda: tk.Label(wind, text='New label').grid(row=3, column=1))
btn3.grid(row=3, column=0)

count = 0
btn4 = tk.Button(wind, text=f'Counter {count}', font=('Times Roman', 16),
                 command=click_counter,
                 activebackground='blue')
btn4.grid(row=4, column=0)

wind.mainloop()