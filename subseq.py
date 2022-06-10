import random
import math
from random import randint

class Subseq:
    def __init__(self, x):
        if isinstance(x, int):
            x = bin(x)[2:]
        self.X = x
    
    def method_1(self, m, k):
        length_of_bits = len(self.X)
        subseq = ''
        for j in range(length_of_bits):
            if j*m+k < length_of_bits:
                subseq += self.X[j*m+k]
        return subseq

    def method_2(self, m):
        length_of_bits = len(self.X)
        num_of_subseq = math.ceil(length_of_bits/m)
        tail = 0
        head = m - 1
        subseq = ''
        for i in range(num_of_subseq):
            if length_of_bits%m == 0:
                s = self.X[tail+i*m:head+i*m]
                subseq += s[randint(0, len(s)-1)]
            else:
                if i == num_of_subseq-1:
                    s = self.X[tail+i*m:head+i*m]
                    subseq += s[randint(0, len(s)-1)]
                else:
                    s = self.X[tail+i*m:len(self.X)-1]
                    subseq += s[randint(0, len(s)-1)]
        return subseq
    
    def method_3(self, p):
        values = [1, 0]
        alpha = random.choices(values, weights=[p, 1-p], k=len(self.X))
        s = ''
        for i in range(len(self.X)):
            if alpha[i] == 1:
                s+=self.X[i]
        return s


    def method_4(self, k):
        values = [1, 0]
        s = ''
        alpha = random.choices(values, weights=[k / len(self.X), (len(self.X)-k)/len(self.X)], k=1)
        for i in range(len(self.X)):
            prob = (k - sum(alpha)) / (len(self.X) - i)
            alpha.append(random.choices(values, weights=[prob, 1 - prob], k=1)[0])

        for i in range(len(self.X)):
            if alpha[i] == 1:
                s+=self.X[i]
        return s

    def method_5(self, m):
        s = ''
        for j in range(len(self.X)-m):
            s+=self.X[j]
        return s