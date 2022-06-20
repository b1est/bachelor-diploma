import math
from scipy.special import erfc as erfc
from scipy.special import gammaincc as gammaincc
from numpy import zeros as zr
from numpy import abs as abs
from numpy import array as array
from numpy import floor as floor
from numpy import max as maxx
from numpy import sqrt as sqrt
from numpy import sum as sum
from scipy.stats import norm as norm

class Tests:
    def __init__(self,  _seq, alpha = 0.01, n = None):
        if isinstance(_seq, str):
            n = len(_seq)
            _seq = int(_seq, 2)
        self.n = n   
        self.seq = _seq
        self.alpha = alpha

    def frequency_test(self):
        s = 0
        for i in range(0, self.n):
            bit = (self.seq >> i) & 1
            s += 2*bit-1
        s_obs = abs(s)/math.sqrt(self.n) 
        p = erfc(s_obs/math.sqrt(2))
        return p>=self.alpha
  
    def frequency_block(self, m = 128):
        N = math.floor(self.n / m)
        if N <= 1:
            return self.frequency_test()
        x_obs = 0
        for i in range(1, N+1):
            pi = 0
            for j in range(1, m+1):
                pi += (self.seq >> ((i-1)*m+j)) & 1 
            pi /= m
            x_obs += math.pow(pi - 0.5, 2)
        x_obs = 4*m*x_obs
        p = gammaincc(N/2, x_obs/2)
        return p>=self.alpha
    
    def runs_test(self):
        pi = 0
        for i in range(self.n):
            pi += (self.seq >> i) & 1
        pi /= self.n
        tau = 2 / math.sqrt(self.n)
        p = 0
        if abs(pi - 0.5) < tau:
            v = 0
            for k in range(self.n-1):
                if (((self.seq & (1 << k)) >> k) & 1 != ((self.seq & (1 << k + 1)) >> k+1) & 1):
                    v += 1
            v += 1
            p = erfc(abs(v-(2*self.n*pi*(1-pi)))/(2*math.sqrt(2*self.n)*pi*(1-pi)))
        return p>=self.alpha

    def longest_one_block_test(self):
        seq = bin(self.seq)[2:]
        while len(seq) < self.n:
            seq = '0' + seq
        if self.n < 128:
            return False
        elif self.n < 6272:
            k = 3
            m = 8
            v = (1, 2, 3, 4)
            pi = (0.2148, 0.3672, 0.2305, 0.1875)
        elif self.n < 750000:
            k = 5
            m = 128
            v = (4, 5, 6, 7, 8, 9)
            pi = (0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124)
        else:
            k = 6
            m = 10000
            v = (10, 11, 12, 13, 14, 15, 16)
            pi = (0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727)
        blocks_number  = math.floor(self.n / m)
        start = 0
        end = m
        x_obs = 0
        freq = zr(k + 1)
        for count in range(blocks_number):
            block_data = seq[start:end]
            max_ones_block = 0
            ones_block_count = 0
            for bit in block_data:
                if bit == '1':
                    ones_block_count += 1
                    max_ones_block = max(max_ones_block, ones_block_count)
                else:
                    max_ones_block = max(max_ones_block, ones_block_count)
                    ones_block_count = 0
            max(max_ones_block, ones_block_count)
            if max_ones_block < v[0]:
                freq[0] += 1
            for j in range(k):
                if max_ones_block == v[j]:
                    freq[j] += 1
            if max_ones_block > v[k - 1]:
                freq[k] += 1
            start += m
            end += m
        for count in range(len(freq)):
            x_obs += pow((freq[count] - (blocks_number  * pi[count])), 2.0) / (blocks_number  * pi[count])
        p = gammaincc(float(k / 2), float(x_obs / 2))
        return p >= self.alpha

    def cumsum_test(self, mode=0):
        seq = bin(self.seq)[2:]
        while len(seq) < self.n:
            seq = '0' + seq
        count = zr(self.n)
        if not mode == 0:
            seq = seq[::-1]
        _count = 0
        for ch in seq:
            sub = 1
            if ch == '0':
                sub = -1
            if _count > 0:
                count[_count] = count[_count -1] + sub
            else:
                count[_count] = sub
            _count += 1
        _abs = maxx(abs(count))
        start = int(floor(0.25 * floor(-self.n / _abs) + 1))
        end = int(floor(0.25 * floor(self.n / _abs) - 1))
        t1 = []
        for k in range(start, end + 1):
            sub = norm.cdf((4 * k - 1) * _abs / sqrt(self.n))
            t1.append(norm.cdf((4 * k + 1) * _abs / sqrt(self.n)) - sub)
        start = int(floor(0.25 * floor(-self.n / _abs - 3)))
        end = int(floor(0.25 * floor(self.n / _abs) - 1))
        t2 = []
        for k in range(start, end + 1):
            sub = norm.cdf((4 * k + 1) * _abs / sqrt(self.n))
            t2.append(norm.cdf((4 * k + 3) * _abs / sqrt(self.n)) - sub)
        p = 1.0 - sum(array(t1))
        p += sum(array(t2))
        return p >= self.alpha
        
    def approximate_entropy_test(self,  block_len=10):
        seq = bin(self.seq)[2:]
        while len(seq) < self.n:
            seq = '0' + seq
        seq += seq[:block_len + 1:]
        max_block = ''
        for i in range(block_len + 2):
            max_block += '1'
        v1_obs = zr(int(max_block[0:block_len:], 2) + 1)
        v2_obs = zr(int(max_block[0:block_len + 1:], 2) + 1)
        for i in range(self.n):
            v1_obs[int(seq[i:i + block_len:], 2)] += 1
            v2_obs[int(seq[i:i + block_len + 1:], 2)] += 1
        v_obs = [v1_obs, v2_obs]
        sums = zr(2)
        for i in range(2):
            for j in range(len(v_obs[i])):
                if v_obs[i][j] > 0:
                    sums[i] += v_obs[i][j] * math.log(v_obs[i][j] / self.n)
        sums /= self.n
        ape = sums[0] - sums[1]
        x_obs = 2.0 * self.n * (math.log(2) - ape)
        p = gammaincc(pow(2, block_len - 1), x_obs / 2.0)
        return p >= self.alpha