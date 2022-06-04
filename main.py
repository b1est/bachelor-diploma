import csv
import random
import numpy as np
from nist import Tests
from subseq import Subseq

print("HEllO")
def embedded_gen(n):
    return "".join(map(str,  np.random.randint(0,2,size=n)))

def Librarian(n): 
    with open("new_text.txt", 'r') as text:
        str_text = text.read()
        s_arr = bytearray(str_text, "utf-8")
        res = ''
        r = random.randint(0, len(str_text) - n)
        for i in range(r, r+n):
            res += bin(s_arr[i])[2:] 
    return res[0:n-1]

def generate_seqs(bits_length, num):
    seq = []
    for i in range(num//2):
        s = embedded_gen(bits_length)
        seq.append(s)
    for i in range(num//2-4):
        s = Librarian(bits_length)
        seq.append(s)
    s = '1'*bits_length
    seq.append(s)
    s = '0'*bits_length
    seq.append(s)
    s = ''
    for i in range(bits_length):
        if i % 2 == 0:
            s += '0'
        else:
            s += '1'
    seq.append(s)
    s = ''
    for i in range(bits_length):
        if i % 2 == 0:
            s += '1'
        else:
            s += '0'
    seq.append(s)
    return seq

def upgrade(s):
    new_s = ''
    for i in range(len(s)):
        if i < len(s)-1:
            new_s += str(int(s[i])^int(s[i+1]))
        else:
            new_s += str(int(s[0])^int(s[len(s)-1]))
    return new_s

def methods_fill(seqs, alpha = 0.05):
    res = []
    for s in seqs:
        ress = []
        for i in range(1, 7):
            if i == 1:
                ress.append(Tests(Subseq(s).method_1(2, 0), alpha).frequency())
            elif i == 2:
                ress.append(Tests(Subseq(s).method_2(4), alpha).frequency())
            elif i == 3:
                ress.append(Tests(Subseq(s).method_3(0.5), alpha).frequency())
            elif i == 4:
                ress.append(Tests(Subseq(s).method_4(len(s)//2), alpha).frequency())
            elif i == 5:
                ress.append(Tests(Subseq(s).method_5(len(s)//2), alpha).frequency())
            else:
                r = nist_tests(s, alpha)
                for re in r:
                    ress.append(re)
        res.append(ress)
    return res

def nist_tests(s, alpha = 0.05):
    res = []
    t = Tests(s, alpha)
    res.append(t.frequency())
    res.append(t.frequency_block())
    res.append(t.runs_test())
    res.append(t.longest_one_block_test())
    res.append(t.cumsum_test())
    res.append(t.cumsum_test(1))
    return res

def make_csv(length, num, name):
    seq = generate_seqs(length, num)
    headers = ["Method 1", "Method 2", "Method 3", "Method 4", "Method 5", "Frequency (Monobit) Test", "Frequency Test within a Block", "Runs Test", "Test for the Longest Run of Ones in a Block", "Cumulative Sums Test (FORWARD)", "Cumulative Sums Test (BACKWARD)"]
    resume = methods_fill(seq)
    with open(name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for r in resume:
            writer.writerow(r)
    return resume

def report_txt(result_tests, name_report_file):
    length_of_seq_result = 0
    headers = [['Kolmogorov True', 'Kolmogorov False', 'NIST True', 'NIST False'], ['GOOD Seq Kolmogorov True', 'GOOD Seq Kolmogorov False', 'GOOD Seq NIST True', 'GOOD Seq NIST False'], ['BAD Seq Kolmogorov True', 'BAD Seq Kolmogorov False', 'BAD Seq NIST True', 'BAD Seq NIST False']]
    with open(name_report_file, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers[0])
        for r in result_tests:
            length_of_seq_result = len(r)
            for i in r:
                t_seq_kolmogorov_count = 0
                f_seq_kolmogorov_count = 0
                t_seq_nist_count = 0
                f_seq_nist_count = 0
                for i in range(length_of_seq_result):
                    if i < 5:
                        if r[i] == True:
                            t_seq_kolmogorov_count += 1
                        if r[i] == False:
                            f_seq_kolmogorov_count += 1
                    else:
                        if r[i] == True:
                            t_seq_nist_count += 1
                        if r[i] == False:
                            f_seq_nist_count += 1
            writer.writerow([t_seq_kolmogorov_count, f_seq_kolmogorov_count, t_seq_nist_count, f_seq_nist_count])
        length_of_result = len(result_tests)
        half = length_of_result//2
        good_seq_tests = [x for l in result_tests[:half] for x in l]
        bad_seq_tests = [x for l in result_tests[half:] for x in l]
        
        def gb_counter(seq):
            t_seq_kolmogorov_count = 0
            f_seq_kolmogorov_count = 0
            t_seq_nist_count = 0
            f_seq_nist_count = 0
            j = 0
            for i in range(len(seq)):
                if i < 5+length_of_seq_result*j:
                    if seq[i] == True:
                        t_seq_kolmogorov_count += 1
                    if seq[i] == False:
                        f_seq_kolmogorov_count += 1
                elif i >= 5+length_of_seq_result*j and i < 11+length_of_seq_result*j:
                    if seq[i] == True:
                        t_seq_nist_count += 1
                    if seq[i] == False:
                        f_seq_nist_count += 1
                    if i == 10+length_of_seq_result*j:
                        j+=1
            
            return [t_seq_kolmogorov_count, f_seq_kolmogorov_count, t_seq_nist_count, f_seq_nist_count]
        
        writer.writerow('\n')
        writer.writerow(headers[1])
        writer.writerow(gb_counter(good_seq_tests))
        writer.writerow('\n')
        writer.writerow(headers[2])
        writer.writerow(gb_counter(bad_seq_tests))

def main():
    num_of_seq = 100
    res256 = make_csv(256, num_of_seq, 'seq256.csv')
    res512 = make_csv(512, num_of_seq, 'seq512.csv')
    report_txt(res256, '256stat.csv')
    report_txt(res512, '512stat.csv')

if __name__ == "__main__": 
    main()