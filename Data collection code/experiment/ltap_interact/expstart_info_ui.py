import tkinter as tkr
class ExpStart(tkr.Frame):
    def __init__(self, master=None):
        tkr.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        tkr.Label(self, text='Now the real experiment will start', font=("Segoe UI", 15, "bold")).pack(ipadx=15, ipady=15)
        tkr.Label(self, text="You can initialize a pause any time the simulation is frozen.", font=("Segoe UI", 12)).pack()
        tkr.Label(self, text="Press 'X' on the keyboard to get in the car.",
                  font=("Segoe UI", 12)).pack()
        tkr.Label(self, text='  ', font=("Segoe UI", 10, "bold")).pack(ipadx=15, ipady=15)
        tkr.Button(self, text='Ok', command=self.proceed, font=("Segoe UI", 10), bg="white").pack()

    def proceed(self):
        self.quit()
