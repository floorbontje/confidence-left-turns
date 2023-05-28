import tkinter as tkr
import random
from tkinter import *



class EndInfoUI(tkr.Frame):
    def __init__(self, master=None):
        tkr.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        tkr.Label(self, text='Thank you for your participation!', font=("Segoe UI", 17, "bold")).pack(ipadx=15, ipady=15)

        tkr.Label(self, text="We only would like to ask you to fill in the additional questionnaire.", font=("Segoe UI", 12)).pack()

        tkr.Label(self, text='  ', font=("Segoe UI", 10, "bold")).pack(ipadx=15, ipady=15)

        tkr.Button(self, text='Ok', command=self.proceed, font=("Segoe UI", 10), bg="white").pack()



    def proceed(self):
        self.quit()

