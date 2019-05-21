import math
import numpy as np
import matplotlib.pyplot as plt


class Element:
    """
    Einzelelement fuer einen Digitalwellenfilter
    Typ 1 fuer 0 > gamma > -0.5
    Typ 2 fuer -0.5 > gamma > -1
    """
    def __init__(self, a):
        self.a = a
        self.delay = np.zeros(2)
        self.ein = 0
        self.aus = 0
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0

    def update(self):
        self.s1 = self.delay[1] - self.ein
        self.s2 = self.a * self.s1 + self.delay[1]
        self.s3 = self.s2 - self.s1

        self.aus = self.s2

    def advance(self):
        self.delay[1] = self.delay[0]
        self.delay[0] = self.s3


class Filter():
    """
    Digitalwellenfilter 19. Grades mit Grenzfrequenz fg.
    Ausgang 1: Tiefpass
    Ausgang 2: Hochpass
    """
    def __init__(self, fs = 16.3e3, F = 64e3):
        self.e = [Element(0)]
        gamma = self.calculate_gamma(fs, F)
        self.e.append(Element(gamma[2]))
        self.e.append(Element(gamma[4]))
        self.e.append(Element(gamma[6]))
        self.e.append(Element(gamma[8]))

        self.e.append(Element(gamma[1]))
        self.e.append(Element(gamma[3]))
        self.e.append(Element(gamma[5]))
        self.e.append(Element(gamma[7]))
        self.e.append(Element(gamma[9]))

        self.aus_hp = 0
        self.auf_tp = 0
        self.ein = 0
        self.delay = 0

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

        self.aus_hp = 0.5 * (self.delay - self.e[9].aus)
        self.aus_tp = 0.5 * (self.delay + self.e[9].aus)

    def advance(self):
        self.delay = self.e[4].aus
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


def test_Filter():
    n = 2**14
    x = np.zeros(n)
    y_tp = np.zeros(n)
    y_hp = np.zeros(n)
    F = 64e3
    fs = 16.3e3
    delta_f = F/n
    fscale = np.arange(0, F/2, delta_f)
    x[0] = 1

    f = Filter(fs, F)

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


test_Filter()
