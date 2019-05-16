import numpy as np
import matplotlib.pyplot as plt


class Element:
    """
    Einzelelement fuer einen Digitalwellenfilter
    Typ 1 fuer 0 > gamma > -0.5
    Typ 2 fuer -0.5 > gamma > -1
    """
    def __init__(self, a, typ=1):
        self.a = a
        self.delay = np.zeros(2)
        self.ein = 0
        self.aus = 0
        self.s1 = 0
        self.s2 = 0
        self.s3 = 0
        if typ in {1, 2, 3}:
            self.typ = typ
        else:
            raise Exception("Ungueltiger Typ '{}'".format(self.typ))

    def update(self):
        if self.typ == 1:
            self.s1 = self.ein + self.delay[1]
            self.s2 = self.s1 * self.a - self.delay[1]
            self.s3 = self.s1 + self.s2

            self.aus = self.s2

        elif self.typ == 2:
            self.s1 = self.ein + self.delay[1]
            self.s2 = self.s1 * self.a - self.delay[1]
            self.s3 = self.s2 - self.s1

            self.aus = self.s3

        elif self.typ == 3:
            self.s1 = self.ein + self.delay[1]
            self.s2 = self.s1 * self.a - self.delay[1]
            self.s3 = self.s1 - self.s2

            self.aus = self.s3

        else:
            raise Exception("Ungueltiger Typ '{}'".format(self.typ))

    def advance(self):
        if self.typ == 1:
            self.delay[1] = self.delay[0]
            self.delay[0] = self.s3

        elif self.typ == 2:
            self.delay[1] = self.delay[0]
            self.delay[0] = self.s2
        elif self.typ == 3:
            self.delay[1] = self.delay[0]
            self.delay[0] = self.s2
        else:
            raise Exception("Ungueltiger Typ '{}'".format(self.typ))


class Filter():
    """
    Digitalwellenfilter 19. Grades mit Grenzfrequenz fg.
    Ausgang 1: Tiefpass
    Ausgang 2: Hochpass
    """
    def __init__(self):
        self.e = [Element(0, 1)]

        self.e.append(Element(0.226119, 1))
        self.e.append(Element(0.397578, 2))
        self.e.append(Element(0.160677, 2))
        self.e.append(Element(0.049153, 2))

        self.e.append(Element(0.063978, 1))
        self.e.append(Element(0.423068, 1))
        self.e.append(Element(0.258673, 2))
        self.e.append(Element(0.094433, 2))
        self.e.append(Element(0.015279, 3))

        self.aus = 0
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

        self.aus = 0.5 * (self.delay + self.e[9].aus)

    def advance(self):
        self.delay = self.e[4].aus
        for e in self.e:
            e.advance()

        # # Ausgang:
        # self.aus = 0.5 * (self.delay - self.e[9].aus)
        # # Oberer Zweig:
        # self.delay = self.e[4].aus
        # self.e[4].ein = self.e[3].aus
        # self.e[3].ein = self.e[2].aus
        # self.e[2].ein = self.e[1].aus
        # self.e[1].ein = self.ein
        # # Unterer Zweig:
        # self.e[9].ein = self.e[8].aus
        # self.e[8].ein = self.e[7].aus
        # self.e[7].ein = self.e[6].aus
        # self.e[6].ein = self.e[5].aus
        # self.e[5].ein = self.ein
        # # Advance the Elements
        # for e_ in self.e:
        #     e_.advance()


def test_Filter():
    n = 2**14
    x = np.zeros(n)
    y = np.zeros(n)
    Fs = 64e3
    delta_f = Fs/n
    fscale = np.arange(0, Fs/2, delta_f)
    x[0] = 1

    f = Filter()

    for index in range(n):
        f.ein = x[index]
        f.update()
        f.advance()
        y[index] = f.aus

    plt.figure(1).set_size_inches(12, 8)
    plt.title("Digitalwellenfilter Zeitbereich")
    plt.plot(x)
    plt.plot(y)
    plt.figure(2).set_size_inches(12, 8)
    plt.title("Digitalwellenfilter Frequenzbereich")
    plt.plot(fscale, -20 * np.log10(abs(np.fft.fft(y)[0:n//2])))
    plt.show()


def test_fft():
    Fs = 64e3
    f1 = 16e3
    f2 = 3e3
    f3 = 31e3
    x = np.zeros(2048)
    y = np.array(range(2048))
    x = (np.sin(2 * np.pi * y * f1 / Fs) + np.sin(2 * np.pi * y * f2 / Fs) +
         np.sin(2 * np.pi * y * f3 / Fs))

    plt.figure(1)
    plt.plot(x)
    plt.figure(2)
    plt.plot(abs(np.fft.fft(x)))
    plt.show()


test_Filter()
