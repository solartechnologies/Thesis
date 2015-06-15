'''
Created on 05/03/2015

@author: Raymond
'''

import serial
from struct import unpack
from analysispack import analysisTools
import numpy as np
from time import sleep
import Tkinter as tk
import ttk
import matplotlib
from Tkconstants import TOP
from math import ceil
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

class uartmicro(serial.Serial):
    
    sampletime = 0
    
    def __init__(self, *args, **kwargs):
        serial.Serial.__init__(self, *args, **kwargs)
        self.baudrate = 250000
        self.timeout = 1
        self.writeTimeout = None
        self.port = 'COM12'
        self.open()
        sleep(1)
        try:
            self.writePIC('h')
            sleep(0.1)
            if(self.read()=='h'):
                print "Connected!"
                self.sampletime = self.sampleTimeADC()
                print self.sampletime*self.getSampleSize()
            else:
                print "Failed connection!"
                self.close()
        except:
            raise
    
    def characterise(self):
        self.writePIC('d')
        self.writePIC('C')
        return self.readADC()
        
    def measureSR(self):
        self.writePIC('d')
        return self.readADC()
        
    def readADC(self):
        """ Queries the microcontroller for data
        and if successful will return a full list, 
        else Exception is returned"""
        
        data = []
        append = data.append
        try:
            wait = self.inWaiting
            read = self.read
            while(wait==0): pass
            while(1):
                values = read(2)
                length = len(values)
                if(length==2):
                    append(unpack('H',values)[0])
                else:
                    if(length==0):
                        break
                    else:
                        raise Exception("Corrupted read: discard data")
                    
                    
            self.flushInput()
        except:
            raise
        
        data = [float(value)*(3.3/1023.0) for value in data]
        return data[:32768]
    
    def sampleTimeADC(self):
        """
        Returns the sample time of the uC ADC module in seconds
        """
        
        pb_t = 25e-9            #peripheral clock time is 25 ns i.e sys_clk/2
        if(not(self.isOpen())):
            return
        
        try:
            self.writePIC('a')
            while(self.inWaiting()==0): pass
            samc = unpack('B', self.read(1))[0]
            adcs = unpack('B', self.read(1))[0]
            self.flushInput()
        except serial.SerialException:
            return serial.SerialException
            
        ad_t = 2*(adcs+1)*pb_t
        aq_t = ad_t*samc
        adcSampleTime = 13*ad_t + aq_t #THE DATASHEET LIES!!!! It takes 13 ADC clock cycles to convert!!
        return adcSampleTime
    
    def getSampleSize(self):
        maxValueRead = 32769
        return maxValueRead
    
    def change_freq(self, frequency, LED):
        """Attempts to change LED frequency"""
        samplePeriod = self.sampletime*self.getSampleSize()

        # sets the frequency to a bin value of the FFT
        if frequency != 0:
            frequency = ceil(frequency*samplePeriod)*(1/float(samplePeriod))
        output = 'f'+ str(frequency) +'f' + str(LED) + 'f' 
        try:
            self.writePIC(output)
            while(self.read()!='f'): pass
        except:
            pass
        return frequency
    
    def writePIC(self, output):
        # easily changed parser to the microcontroller
        try:
            self.write(output)
            self.write('\n')
        except:
            pass


class startPage(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        self.title('Raytech')
        self.minsize(width=600, height=400)
        
        optionbar = self.optionbar_init()
        self.config(menu=optionbar)
        
        # Set up frames for values and options
        container = tk.Frame(self)
        container.pack(side='top', fill=tk.BOTH, expand=True, pady=2, padx=2)        
        
        optionsLabel = ttk.Label(container, text="Options", font=("purisa", 16))
        options = ttk.LabelFrame(container, labelwidget = optionsLabel)
        settings = ttk.LabelFrame(container, height = 50, text='Settings')
        self.option_frame_init(options)
        
        pages = ttk.Notebook(container)
        self.note = pages
        plusTab = tk.Frame(pages)
        initGraph = graphTab(pages)
        self.graphbox = [initGraph]
        initGraph.pack(fill='both', expand=True)
        
        pages.add(initGraph, text="  Data  ", padding=1)
        pages.add(plusTab,text='+')
        
        pages.bind("<<NotebookTabChanged>>", lambda e: self.addTab())
        
        """ configuring the page grid """
        options.grid(row=1, column=0, sticky ='nesw')
        settings.grid(row=2, column=0, columnspan=2, sticky='nsew')
        pages.grid(row=1,column=1, sticky='nsew')
        container.grid_columnconfigure(0, weight=1, minsize=60)
        container.grid_columnconfigure(1, weight=2)
        container.grid_rowconfigure(0, minsize = 1)
        container.grid_rowconfigure(1, weight=2)
        container.grid_rowconfigure(2, weight=1)
    
    def option_frame_init(self, options):
        # buttons to put in options frame
        button_new = ttk.Button(options, text="New Measurement", 
                                        command = self.measureCell)
        button_load = ttk.Button(options, text="Load EQE")
        button_save = ttk.Button(options, text="Save Measurement")
        
        button_new.pack(padx = [2,5], pady = [8,2], expand=False, fill=tk.X, side = TOP)
        button_load.pack(padx = [2, 5], pady = 2, expand=False, fill=tk.X, side = TOP)
        button_save.pack(padx = [2, 5], pady = 2, expand=False, fill=tk.X, side = TOP)
        
    def optionbar_init(self):
        """
        Creates the drop down file menu for the GUI
        """
        optionbar = tk.Menu(self)
        
        filemenu = tk.Menu(optionbar, tearoff=0)
        filemenu.add_command(label='Open')
        filemenu.add_command(label='Save')
        filemenu.add_command(label='Save_as')
        filemenu.add_separator()
        filemenu.add_command(label='Exit')
        optionbar.add_cascade(label="File", menu=filemenu)
        
        configmenu = tk.Menu(optionbar, tearoff=0)
        configmenu.add_command(label="Run Diagnostics")
        configmenu.add_command(label="Settings")
        optionbar.add_cascade(label='Config', menu=configmenu)
        
        helpmenu = tk.Menu(optionbar, tearoff=0)
        helpmenu.add_command(label="About")
        optionbar.add_cascade(label="Help", menu=helpmenu)
        
        return optionbar
    
    def addTab(self):
        slot = self.note.index('end')-1
        if slot == self.note.index(self.note.select()):
            tab = graphTab(self.note)
            self.toolbox.append(tab)
            tab.pack(fill='both', expand=True)
            self.note.insert(slot, tab,text="  Data {0} ".format(slot))
            self.note.select(slot)
    
    def plotData(self, x, y):
        currTab = self.getTab()
        currTab.plotData(x,y)
    
    def measureCell(self):
        # TODO make this vialble for multiple LED's
        # TODO move this to graphTab and just call the method there
        currTab = self.getTab()
        stream = uartmicro()
        stream.change_freq(currTab.analysis.LED['585nm']['frequency'])
        rawdata = stream.readADC()
        currTab.analysis.rawData = rawdata
        fx, fy = currTab.analysis.fastFourierDetrend(stream.sampleTimeADC()*1e9)
        currTab.analysis.amplitudeLED()
        
        
    def getTab(self):
        currTab = self.graphbox[self.note.index("current")]
        return currTab
        
    def exit(self):
        pass
    
class sr_settings(ttk.Frame):

    def __init__(self, controller):
        pass

class graphTab(tk.Frame):
    
    # TODO pass a copy of
    analysis = analysisTools()
    def __init__(self, parent):
        tk.Frame.__init__(self)
        
        f = Figure(figsize = (5,4),dpi=100)
        self.graph = f.add_subplot(111)
        
        canvas = FigureCanvasTkAgg(f,self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM,fill=tk.BOTH,expand=True)
        
        toolbar = NavigationToolbar2TkAgg(canvas,self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP,fill=tk.BOTH, expand= True)
    
    def plotData(self, x, y):
        self.graph.plot(x,y)




    

        
        

windowOne = startPage()
windowOne.mainloop()

        
"""
def update():
    
    windowOne.after(1000,update)
"""
# data = []
# with open('test.csv','r') as f:
#     iterate = f
#     iterate.readline()
#     line = iterate.readline()
#     line = line.split(',')
#     sampletime = float(line[0])
#     for line in f:
#         line = line.split(',')
#         data.append(float(line[1]))
"""
stream = uartmicro()
analysis = analysisTools()
for freq, LED in analysis.getLEDfreq():
    print stream.change_freq(freq, LED)
# data = stream.measureSR()
# data_len = len(data)
# print data_len
sampletime = stream.sampleTimeADC()
data_ref = stream.characterise()
analysis.setData(data_ref, sampletime)
analysis.characterise()

fx, fy = analysis.fastFourier(sampletime, data_ref)
fft_len = len(fy)
with open('test1.csv','w') as f:
    for index in xrange(len(data_ref)-1):
        f.write(str(index*sampletime)+',')
        f.write(str(+data_ref[index])+",")
        if(index < fft_len):
            f.write(str(fx[index])+","+
                    str(fy[index])+",")
        f.write("\n")


raw_input('ready?')

for freq, LED in analysis.getLEDfreq():
    print stream.change_freq(freq, LED)
data_test = stream.measureSR()
analysis.setData(data_test, sampletime)
analysis.quantumEfficiency()
wave, EQE = analysis.getEQE()

with open('result.csv', 'w') as f:
    
    for idx in range(len(wave)):
        f.write(str(wave[idx])+',' + str(EQE[idx])+'\n')
        print str(wave[idx]) + ': ' + str(EQE[idx])

fx, fy = analysis.fastFourier(sampletime, data_test)
fft_len = len(fy)
with open('test2.csv','w') as f:
    for index in xrange(len(data_test)-1):
        f.write(str(index*sampletime)+',')
        f.write(str(+data_test[index])+",")
        if(index < fft_len):
            f.write(str(fx[index])+","+
                    str(fy[index])+",")
        f.write("\n")

# fx, fy = analysis.fastFourier(sampletime, data)
# fft_len = len(fy)
# print fft_len
# with open('test.csv','w') as f:
#     for index in xrange(0,len(data)-1):
#         f.write(str(index*sampletime)+',')
#         f.write(str(+data[index])+",")
#         if(index < fft_len):
#             f.write(str(fx[index])+","+
#                     str(fy[index])+",")
#         f.write("\n")
#     f.close()
stream.close()

"""
"""
windowOne = EQE_TkGui.startPage()
windowOne.after(1000, update)
windowOne.plotData( fx, fy)
windowOne.mainloop()
"""
