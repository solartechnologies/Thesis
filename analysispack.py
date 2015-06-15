import numpy as np
from scipy import signal, integrate
from scipy.interpolate import interp1d
from math import ceil, log

class analysisTools:
    rawData = {'time':[],'voltage':[]}        # holds the data read from the ADC
    
    LEDreference = {}                       # This is where all reference and 
    EQEreff = {}                            # characterisation values are held
    planks_constant = 6.62606957e-34
    speed_of_light = 299792458
    electron_charge = 1.602177e-19
    
    EQE = {'EQE':[],'Wavelength':[]}         # Measured value for the cell
    
    def __init__(self):
        if not self.EQEreff:
            self.loadRefCell()
        
        if not self.LEDreference:
            self.setLEDReference()
    
    def setData(self, data, sampletime):
        self.rawData['time'] = [sampletime*num for num in range(len(data))]
        self.rawData['voltage'] = data
        
    def getLEDfreq(self):
        """creates an iterable function that returns duel
        values frequency and Channel"""
        for key, content in self.LEDreference.iteritems():
            yield [content['frequency'], content['Channel']]
    
    def setLEDReference(self):
        # stores the frequency, integrals constants and scalars
        # for each LED in file
        freq = 6000             # start at 10 kHz for stable gain from the TIA
        LEDdata = self.loadLEDspectra()
        LEDref = {}
        trapz = integrate.trapz
        for key, content in LEDdata.iteritems():
            LEDref[key] = {'frequency': freq}
            update = LEDref[key].update
            update(content)
            freq += 580
            
            integral = trapz(content['SpectraCounts'], content['Wavelength'])
            update({'Integral':integral})
            update({'Scalar': 1.0/integral})
            norm = []
            for idx in xrange(len(content['Wavelength'])):
                norm.append(content['SpectraCounts'][idx]*content['Wavelength'][idx]*1e-9)
            integral = trapz(norm, content['Wavelength'])
            update({'spectral_integral': integral})
            
            interp = interp1d(self.EQEreff['Wavelength'], self.EQEreff['EQE'])
            min_wave = min(self.EQEreff['Wavelength'])
            max_wave = max(self.EQEreff['Wavelength'])
            Wavelength = []
            append_Wave = Wavelength.append
            ratio = []
            append_ratio = ratio.append
            for idx, omega in enumerate(content['Wavelength']):
                if omega < max_wave and omega > min_wave:
                    append_Wave(omega)
                    append_ratio(interp(omega)*content['SpectraCounts'][idx]*omega*1e-9)
            integral = trapz(ratio, Wavelength)
            update({'spectral_eqe_integral': integral})
            update(content)
        
        self.LEDreference.update(LEDref) #update original class variable to be shared across instances
    
    def getEQE(self):
        # returns the wavelength and measured eqe
        return self.EQE['Wavelength'], self.EQE['EQE']*100
    
    def quantumEfficiency(self):
        current = self.currentMeasure(2.5e5)
        eqe = []
        wavelength = []
        for key in current:
            if not self.LEDreference[key]['Pmax']:
                raise Exception('LEDs must be characterised first')
            
            scalar = self.LEDreference[key]['Scalar']
            numerator = current[key]*self.planks_constant*self.speed_of_light
            denomenator = self.electron_charge*scalar* \
                            self.LEDreference[key]['spectral_integral']* \
                            self.LEDreference[key]['Pmax']
            eqe.append(numerator/denomenator)
            
            # need to use the median wavelength value not the peak
            
            medi = [R for R in self.LEDreference[key]['SpectraCounts'] \
                     if ((R*scalar)>0.01) ]
            length = len(medi)
            if length%2 != 0:
                result = medi[length/2]
            else:
                odd = length/2 -1
                even = length/2
                result = float(medi[even] + medi[odd])/2.0
            print result
            idx = findIndx(self.LEDreference[key]['SpectraCounts'],result)
            wavelength.append(self.LEDreference[key]['Wavelength'][idx])
        
        # sort wavelength and eqe data so that they wavelength is sorted
        zipped = zip(wavelength, eqe)
        zipped.sort()
        self.EQE['Wavelength'], self.EQE['EQE'] = zip(*zipped)
        
    def characterise(self):
        current = self.currentMeasure(2.5e5)
        
        for key in current:
            numerator = current[key]*self.planks_constant*self.speed_of_light
            denomenator = self.electron_charge*self.LEDreference[key]['Scalar']*\
                                self.LEDreference[key]['spectral_eqe_integral']
            pmax = numerator/denomenator
            self.LEDreference[key].update({'Pmax': pmax})
    
    def loadRefCell(self):
        # Loads the EQE distribution from  for the reference cell used
        # in the current measurements
        with open('EQEreff.csv','r') as f:
            
            EQE = {'Wavelength': [],'EQE': []}
            append_Wave = EQE['Wavelength'].append
            append_EQE = EQE['EQE'].append
            next(f)
            for line in f:
                line = line.strip('\n')
                line = line.split(',')
                append_Wave(float(line[0]))
                append_EQE(float(line[1]))
            self.EQEreff.update(EQE)
        
    @staticmethod
    def loadLEDspectra():
        # For each intensity normalised spectral distribution contained in 
        # LEDspectra.csv load the values into a dictionary
        # New spectra can be added by pasting into the file with the same
        # format
        with open('LEDspectra.csv','r') as f:
            
            keys = f.readline()
            keys = keys.strip('\n')
            keys = keys.split(',')
            channel = keys[1::2]
            keys = keys[::2]
            LED = {}
            for idx, key in enumerate(keys):
                LED[key] = {'Wavelength':[],"SpectraCounts": [],\
                            'Channel': channel[idx]} 
            for line in f:
                line = line.strip('\n')
                line = line.split(',')
                for idx, key in enumerate(keys):
                    LED[key]['Wavelength'].append(float(line[idx*2]))
                    LED[key]['SpectraCounts'].append(float(line[idx*2+1]))
            return LED
    
    def fastFourier(self, t_sample, data):
        """
        >>> t_sample is the sample time in seconds
        >>> data is the periodic data to be transformed
        >>> the data is detrended to remove linear (DC) component
        """

        detrended = signal.detrend(data)            #remove DC signal
        length = len(detrended)                    # save data length
        # find the next power of two for efficient FFT calculation
        NFFT = 2**int((ceil(log(length,2)))) 
        fy = np.fft.fft(detrended,NFFT)         #calculate fft
        fy = 2.0*abs(fy)/NFFT               #normalise the real component
        fx = np.fft.fftfreq(NFFT, t_sample)
        return fx[:NFFT/2],fy[:NFFT/2]
    
    def currentMeasure(self, gain):
        """
        Returns the current measured by the ADC
        """
        if not self.rawData['voltage']:
            raise Exception('Invalid Operation: No data has been found')
        
        fftx, ffty = self.fastFourier(self.rawData['time'][1], self.rawData['voltage'])
        LED = self.LEDreference
        current = {}
        update = current.update
        for key in LED:
            idx = findIndx(fftx, LED[key]['frequency'])
            amplitude = 2*ffty[idx]
            update({key:amplitude/float(gain)})
        
        return current
            
            

def findIndx(vals, value):
    least = abs(vals[0] - value)
    idx = None
    for i,num in enumerate(vals):
        if num == value:
            return i
        elif abs(num-value) < least:
            least = abs(num-value)
            idx = i
    return idx