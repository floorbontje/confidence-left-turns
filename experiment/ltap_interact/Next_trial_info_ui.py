import tkinter as tkr
from tkinter import *



class NextInfoUI(tkr.Frame):
    def __init__(self, master=None):
        tkr.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        self.question = tkr.Label(self, text='Ready for the next route?', font=("Segoe UI", 17, "bold")).pack(ipadx=15, ipady=15)

        self.yes_button = tkr.Button(self, text='Yes', command=self.proceed, font=("Segoe UI", 12), bg="white")
        self.yes_button.pack()

        self.white = tkr.Label(self, text='  ', font=("Segoe UI", 10, "bold")).pack()

        self.no_button = tkr.Button(self, text='No', command=self.proceed_no, font=("Segoe UI", 12), bg="white")
        self.no_button.pack()

        tkr.Label(self, text='Press "X" on the keyboard to get in the car', font=("Segoe UI", 12)).pack(ipadx=15, ipady=15)



    def proceed(self):
        self.next_info = {'next': 1}
        self.quit()

    def proceed_no(self):
        self.next_info = {'next': 0}
        self.quit()