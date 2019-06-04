import math
import numpy as np
import matplotlib.pyplot as plt


class Element:
    """
    Einzelelement fuer einen Digitalwellenfilter
    Typ 1 fuer 0 > gamma > -0.5
    Typ 2 fuer -0.5 > gamma > -1
    """
    def __init__(self, a, T=1):
        self.T = T
        self.a = a
        self.delay = np.zeros(2*self.T)
        self.ein = 0
        self.aus = 0
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0

    def update(self):
        self.s1 = self.delay[2*self.T-1] - self.ein
        self.s2 = self.a * self.s1 + self.delay[2*self.T-1]
        self.s3 = self.s2 - self.s1

        self.aus = self.s2

    def advance(self):
        for i in reversed(range(1,2*self.T)):
            self.delay[i] = self.delay[i-1]

        self.delay[0] = self.s3


class Filter():
    """
    Digitalwellenfilter 19. Grades mit Grenzfrequenz fg.
    Ausgang 1: Tiefpass
    Ausgang 2: Hochpass
    """
    def __init__(self, fs=16.3e3, F=64e3, T=1):
        self.T = T
        self.e = [Element(0)]
        gamma = self.calculate_gamma(fs, F)
        self.e.append(Element(gamma[2], self.T))
        self.e.append(Element(gamma[4], self.T))
        self.e.append(Element(gamma[6], self.T))
        self.e.append(Element(gamma[8], self.T))

        self.e.append(Element(gamma[1], self.T))
        self.e.append(Element(gamma[3], self.T))
        self.e.append(Element(gamma[5], self.T))
        self.e.append(Element(gamma[7], self.T))
        self.e.append(Element(gamma[9], self.T))

        self.aus_hp = 0
        self.auf_tp = 0
        self.ein = 0
        self.delay = np.zeros(self.T)

    def update(self):
        self.e[1].ein = self.ein
        self.e[5].ein = self.ein
        self.e[1].update()
        self.e[5].update()

        self.e[2].ein = self.e[1].aus
        self.e[6].ein = self.e[5].aus
        self.e[2].update()
        self.e[6].update()

        self.e[3].ein = self.e[2].aus
        self.e[7].ein = self.e[6].aus
        self.e[3].update()
        self.e[7].update()

        self.e[4].ein = self.e[3].aus
        self.e[8].ein = self.e[7].aus
        self.e[4].update()
        self.e[8].update()

        self.e[9].ein = self.e[8].aus
        self.e[9].update()

        self.aus_hp = 0.5 * (self.delay[self.T-1] - self.e[9].aus)
        self.aus_tp = 0.5 * (self.delay[self.T-1] + self.e[9].aus)

    def advance(self):
        for i in reversed(range(1,self.T)):
            self.delay[i] = self.delay[i-1]
            
        self.delay[0] = self.e[4].aus
        
        for e in self.e:
            e.advance()


    def calculate_gamma(self,fs, F):
        gamma = np.zeros(10)
        q = np.zeros(5)
        q[0] = math.tan(math.pi * fs / F)

        for i in range(1,5):
            q[i] = q[i-1]**2 + (q[i-1]**4 - 1)**(1/2)

        c = np.zeros(5)
        y = np.zeros(10)
        for i in range(1,10):
            c[4] = q[4]/(math.sin(i * math.pi / 19))
            for k in reversed(range(0,4)):
                c[k] = 1 / (2 * q[k]) * (c[k+1] + 1/c[k+1])
                
            y[i] = 1/c[0]

        a = np.zeros(10)

        for i in range(1,10):
            a[i] = 2 / (1 + y[i]**2) * (1 - (q[0]**2 + 1 / (q[0]**2) - y[i]**2) * y[i] **2)**(1/2)
            gamma[i] = (a[i] - 2) / (a[i] + 2)
        
        return gamma

class Filter_Bank():
    """
    Oktavfilterbank mit 9 Ausgaengen. Hoechste Frequenz an aus[0], niedrigste an aus[8]
    """

    def __init__(self, fs = 16.3e3, F = 64e3, delay_num = [1, 1, 1, 1, 1, 1, 1, 1, 1]):
        self.delay = []
        for i in delay_num:
            self.delay.append(np.zeros(i))
        self.F = []      
        for i in range(8):
            self.F.append(Filter(fs, F, 2**i))
        self.ein = 0
        self.aus = np.zeros(9)

    def update(self):
        self.F[0].ein = self.ein
        self.F[0].update()
        self.delay[0][0] = self.F[0].aus_hp

        for i in range(1,7):
            self.F[i].ein = self.F[i-1].aus_tp
            self.F[i].update()
            self.delay[i][0] = self.F[i].aus_hp

        self.F[7].ein = self.F[6].aus_tp
        self.F[7].update()
        self.delay[7][0] = self.F[7].aus_hp

        self.delay[8][0] = self.F[7].aus_tp


    def advance(self):
        for f in self.F:
            f.advance()

        for i in range(9):
            for j in reversed(range(1,len(self.delay[i]))):
                self.delay[i][j] = self.delay[i][j-1]
            self.aus[i] = self.delay[i][-1]
        
        


def test_Filter(T=1):
    n = 2**14
    x = np.zeros(n)
    y_tp = np.zeros(n)
    y_hp = np.zeros(n)
    F = 64e3
    fs = 16.3e3
    delta_f = F/n
    fscale = np.arange(0, F/2, delta_f)
    x[0] = 1

    f = Filter(fs, F, T)

    for index in range(n):
        f.ein = x[index]
        f.update()        
        f.advance()
        y_tp[index] = f.aus_tp
        y_hp[index] = f.aus_hp
    plt.figure(1).set_size_inches(12, 8)
    plt.title("Digitalwellenfilter Frequenzbereich")
    plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(y_tp)[0:n//2])))
    plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(y_hp)[0:n//2])))
    plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(y_hp + y_tp)[0:n//2])))    
    plt.show()

def test_Filter_Bank(delay_num = [1, 1, 1, 1, 1, 1, 1, 1, 1]):
    n = 2**14
    x = np.zeros(n)
    y = np.zeros([9,n])
    F = 64e3 * 2/3
    fs = 16.3e3 * 2/3
    delta_f = F/n
    fscale = np.arange(0, F/2, delta_f)
    x[0] = 1

    b = Filter_Bank(delay_num = delay_num)

    for index in range(n):
        b.ein = x[index]
        b.update()
        b.advance()
        for i in range(9):
            y[i-1][index] = b.aus[i-1]

    # for i in range(9):
    #     plt.figure(i).set_size_inches(12, 8)
    #     plt.title("Digitalwellenfilterbank Daempfung im Frequenzbereich Ausgang {0}".format(i))
    #     plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(y[i])[0:n//2])))
        
    plt.figure(10).set_size_inches(12, 8)
    plt.title("Digitalwellenfilterbank Daempfung im Frequenzbereich Summe der Ausgaenge")
    plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(sum(y))[0:n//2])))

    plt.figure(11).set_size_inches(12, 8)
    plt.title("Digitalwellenfilterbank Daempfung im Zeitbereich Alle Augaenge")
    for i in range(9):
        plt.plot(y[i])

    # plt.legend([str(x) for x in range(9)])    
    plt.show()
    
# test_Filter(1)

maximum = [4, 20, 44, 91, 186, 375, 755, 1500, 900]

delay_num = [max(maximum) - m if max(maximum) - m > 0 else 1 for m in maximum]
test_Filter_Bank(delay_num)
