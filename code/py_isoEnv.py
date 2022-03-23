import numpy as np

def iso_distr(temp):
    hydrogen = int(temp[1])

    carbon = int(temp[0])

    nitrogen = int(temp[2])

    oxygen = int(temp[3])

    sulfur = int(temp[4])

    pH = [0.999885, 0.0001157]
    pC = [0.9893, 0.0107]
    pN = [0.99632, 0.00368]
    pO = [0.99757, 0.00038, 0.00205]
    pS = [0.9493, 0.0076, 0.0429, 0.0002]

    p = convolve(carbon, pC)
    p = np.convolve(p, convolve(oxygen, pO))
    p = np.convolve(p, convolve(hydrogen, pH))
    p = np.convolve(p, convolve(nitrogen, pN))
    p = np.convolve(p, convolve(sulfur, pS))
    
    iso = np.array(cut(p / np.max(p)),dtype="float64")
    return iso


def bits1(n):
    b = []
    while n:
        b = [n & 1] + b
        n >>= 1
    return b or [0]


def convolve(number, probability):
    bitarray = bits1(number)
    pi = probability
    p = [1]
    for i, b in enumerate(bitarray[::-1]):
        p = cut(np.convolve(p, pi)) if b == 1 else p
        pi = cut(np.convolve(pi, pi))

    return p


def cut(array,tr=0.00001):

    index = np.where(array > tr)[0][-1]

    if (len(array) > index):
        return array[:index + 1]
    else:
        return (array)
      
      
