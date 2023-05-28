import tkinter as tkr
import random
from tkinter import *



class StartInfoUI(tkr.Frame):
    def __init__(self, master=None):
        tkr.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        tkr.Label(self, text='Welkom', font=("Segoe UI", 17, "bold")).pack(ipadx=15, ipady=15)

        tkr.Label(self, text="1. Press 'X' on the keyboard to get into the car", font=("Segoe UI", 10)).pack(ipadx=15, ipady=15)
        tkr.Label(self, text="2. Indicate your initially planned driving behaviour each time have to make a left turn", font=("Segoe UI", 10)).pack()
        tkr.Label(self, text="  'Go': O-button "
                             "  'Wait': right handle ", font=("Segoe UI", 10)).pack()

        tkr.Label(self, text='  ', font=("Segoe UI", 10, "bold")).pack(ipadx=5, ipady=5)

        tkr.Button(self, text='Ok', command=self.proceed, font=("Segoe UI", 10), bg="white").pack()



    def proceed(self):
        self.quit()

