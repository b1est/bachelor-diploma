import csv
import random
import numpy as np
from nist import Tests
from subseq import Subseq

import time
import os
import glob
from xlsxwriter.workbook import Workbook

def embedded_gen(n):
    return "".join(map(str, np.random.randint(0,2,size=n)))

def Librarian(n, r): 
    with open("new_text.txt", 'r') as text:
        str_text = text.read()
        s_arr = bytearray(str_text, "utf-8")
        res = ''
       
        for i in range(r, r+n):
            byte = bin(s_arr[i])[2:]
            while len(byte) < 8:
                byte = '0' + byte
            res += byte
    return res

def generate_seqs(bits_length, num):
    seq = []
    for i in range(num//2):
        s = embedded_gen(bits_length)
        seq.append(s)
    byte_length = bits_length//8
    s = open("new_text.txt", 'r')
    text_len = len(s.read())
    s.close()
    R = []
    for i in range(num//2-4):
        r = random.randint(0, text_len - byte_length)

        while r in R and r in [ri for ri in range(r, r+byte_length)]:
            r = random.randint(0, text_len - byte_length)
        s = Librarian(byte_length, r)
        
        for ri in range(r, r+byte_length):
            R.append(ri)
        seq.append(s)
    s = '1'*bits_length
    seq.append(s)
    s = '01'*(bits_length//2)
    seq.append(s)
    s = '0011'*(bits_length//4)
    seq.append(s)
    s = '0000000011111111'*(bits_length//16)
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
        si = Subseq(s)
        s_len = len(s)
        for i in range(1, 7):
            if i == 1:
                r = si.method_1(2, 0)
                ress.append(Tests(r, alpha).frequency())
                r = si.method_1(2, 1)
                ress.append(Tests(r, alpha).frequency())
            elif i == 2:
                if s_len == 256:
                    r = si.method_2(2)
                elif s_len == 512:
                    r = si.method_2(4)
                ress.append(Tests(r, alpha).frequency())
            elif i == 3:
                r = si.method_3(0.5)
                
                ress.append(Tests(r, alpha).frequency())
            elif i == 4:
                r = si.method_4(len(s)//2)
                
                ress.append(Tests(r, alpha).frequency())
            elif i == 5:
                r = si.method_5(len(s)//2)
                
                ress.append(Tests(r, alpha).frequency())
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

def make_csv(resume, name):
    

    headers = ["№", "Method 1a", "Method 1b", "Method 2", "Method 3", "Method 4", "Method 5", "Frequency (Monobit) Test", "Frequency Test within a Block", "Runs Test", "Test for the Longest Run of Ones in a Block", "Cumulative Sums Test (FORWARD)", "Cumulative Sums Test (BACKWARD)"]
    
    with open(name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        i = 1
        for r in resume:
            tmp = [i]
            for ri in r:
                tmp.append(ri)
            writer.writerow(tmp)
            i+=1
    return resume

def report(result_tests, name_report_file):
    length_of_seq_result = 0
    headers = [["№",'Kolmogorov True', 'Kolmogorov False', 'NIST True', 'NIST False'], ['GOOD Seq Kolmogorov True', 'GOOD Seq Kolmogorov False', 'GOOD Seq NIST True', 'GOOD Seq NIST False'], ['BAD Seq Kolmogorov True', 'BAD Seq Kolmogorov False', 'BAD Seq NIST True', 'BAD Seq NIST False']]
    with open(name_report_file, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers[0])
        ii = 1
        for r in result_tests:
            length_of_seq_result = len(r)
            for i in r:
                t_seq_kolmogorov_count = 0
                f_seq_kolmogorov_count = 0
                t_seq_nist_count = 0
                f_seq_nist_count = 0
                for i in range(length_of_seq_result):
                    if i < 6:
                        if r[i] == True:
                            t_seq_kolmogorov_count += 1
                        if r[i] == False:
                            f_seq_kolmogorov_count += 1
                    else:
                        if r[i] == True:
                            t_seq_nist_count += 1
                        if r[i] == False:
                            f_seq_nist_count += 1
            writer.writerow([ii, t_seq_kolmogorov_count, f_seq_kolmogorov_count, t_seq_nist_count, f_seq_nist_count])
            ii+=1
        length_of_result = len(result_tests)
        half = length_of_result//2
        good_seq_tests = [x for l in result_tests[half:] for x in l]
        bad_seq_tests = [x for l in result_tests[:half] for x in l]
        def gb_counter(seq):
            t_seq_kolmogorov_count = 0
            f_seq_kolmogorov_count = 0
            t_seq_nist_count = 0
            f_seq_nist_count = 0
            j = 0
            for i in range(len(seq)):
                if i < 6+length_of_seq_result*j:
                    if seq[i] == True:
                        t_seq_kolmogorov_count += 1
                    if seq[i] == False:
                        f_seq_kolmogorov_count += 1
                elif i >= 6+length_of_seq_result*j and i < 12+length_of_seq_result*j:
                    if seq[i] == True:
                        t_seq_nist_count += 1
                    if seq[i] == False:
                        f_seq_nist_count += 1
                    if i == 11+length_of_seq_result*j:
                        j+=1
            return [t_seq_kolmogorov_count, f_seq_kolmogorov_count, t_seq_nist_count, f_seq_nist_count]    
        writer.writerow('\n')
        writer.writerow(headers[1])
        writer.writerow(gb_counter(bad_seq_tests))
        writer.writerow('\n')
        writer.writerow(headers[2])
        writer.writerow(gb_counter(good_seq_tests))

def convert_csv_to_xlsx():
    for csvfile in glob.glob(os.path.join('.', '*.csv')):
        workbook = Workbook(csvfile[:-4] + '.xlsx')
        worksheet = workbook.add_worksheet()
        with open(csvfile, 'rt', encoding='utf8') as f:
            reader = csv.reader(f)
            for r, row in enumerate(reader):
                for c, col in enumerate(row):
                    worksheet.write(r, c, col)
        workbook.close()

def clear_folder_from_file_type(_type):
    dir_list = os.listdir()
    type_len = -1*len(_type)
    for file in dir_list:
        if file[type_len:] == _type:
            os.remove(file)

def main(bits_len = 256, num_of_seq = 100):
    alpha = [0.01, 0.05, 0.1]
    bits_len = [256, 512]
    for bl in bits_len:
        print(f'Довжина {bl} біт')
        seq = generate_seqs(bl, num_of_seq)
        for a in alpha:
            print(f'\tРівень значущості {a}')
            resume = methods_fill(seq, a)
            res = make_csv(resume, f'Результати тестів послідовності довжиною {bl} біт({int(a*100)}%).csv')
            report(res, f'{bl}.{int(a*100)}.stat.csv')
    
    
if __name__ == "__main__": 
    clear_folder_from_file_type('xlsx')
    main()
    convert_csv_to_xlsx()
    time.sleep(3)
    clear_folder_from_file_type('csv')