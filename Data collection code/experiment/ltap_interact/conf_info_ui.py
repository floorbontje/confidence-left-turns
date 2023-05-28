import tkinter as tkr
from tkinter import *



class ConfInfoUI(tkr.Frame):
    def __init__(self, master=None):
        tkr.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()


    def createWidgets(self):

        tkr.Label(self, text='How confident were you in your decision?',
                  font=("Segoe UI", 15, "bold")).pack(ipadx=15, ipady=15)
        tkr.Label(self, text='Give a rate on a scale of 1 (not at all) to 5 (very)', font=("Segoe UI", 12)).pack(ipadx=15, ipady=15)


        frame = LabelFrame(self)
        frame.pack(pady=30)
        self.conf_var = tkr.IntVar(value=1)
        conf_radios = [tkr.Radiobutton(frame, text=str(rating), padx=20,
                                       variable=self.conf_var, value=rating)
                       for rating in [1, 2, 3, 4, 5]]

        for rating in range(0, 5):
            conf_radios[rating].grid(row=5, column=rating)

        tkr.Label(frame, text='Not so', font=("Segoe UI", 8)).grid(row=4, column=0)
        tkr.Label(frame, text='Very', font=("Segoe UI", 8)).grid(row=4, column=4)


        self.submit_button = tkr.Button(self, text='Submit', command=self.proceed, font=("Segoe UI", 10), bg="white")
        self.submit_button.pack()

        tkr.Label(self, text='To continue driving press the "L3" button on the steering wheel',
                  font=("Segoe UI", 12)).pack(ipadx=15, ipady=15)

    def proceed(self):
        self.conf_info = {'conf': self.conf_var.get()}
        self.quit()
