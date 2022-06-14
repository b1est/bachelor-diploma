import csv
import random
import numpy as np
from nist import Tests
from subseq import Subseq
from confident import confidence_interval_probability_bernoulli_model_independent_examinations as cipbmie

import time
import os
import glob
from xlsxwriter.workbook import Workbook

def embedded_gen(n):
    return "".join(map(str, np.random.randint(0,2,size=n)))

def bad_generator(p1, p2, bits_length, size=1):
    res = []
    if isinstance(p1, int) and isinstance(p2, int):
        for i in range(size):
            s = ''
            seq = random.choices([1, 0], weights=[p1, p2], k=bits_length)
            for si in seq:
                s += str(si)
            
            res.append(s)
    else:
        if size % 2 == 0:
            for j in range(size):
                if j == size//2 - 1:
                    s = ''
                    seq = random.choices([1, 0], weights=[p1[0], p2[0]], k=bits_length)
                    for si in seq:
                        s += str(si)
                    
                    res.append(s)
                else:
                    s = ''
                    seq = random.choices([1, 0], weights=[p1[1], p2[1]], k=bits_length)
                    for si in seq:
                        s += str(si)
                    
                    res.append(s)
        else:
            for j in range(size):
                if j == size//2+1:
                    s = ''
                    seq = random.choices([1, 0], weights=[p1[0], p2[0]], k=bits_length)
                    for si in seq:
                        s += str(si)
                    
                    res.append(s)
                else:
                    s = ''
                    seq = random.choices([1, 0], weights=[p1[1], p2[1]], k=bits_length)
                    for si in seq:
                        s += str(si)
                    
                    res.append(s)
    return res
    
def generate_seqs(bits_length, num):
    seq = []
    for i in range(num//2):
        s = embedded_gen(bits_length)
        seq.append(s)
    s = bad_generator((0.6, 0.8), (0.4, 0.2), bits_length, num//2)
    for item in s:
        seq.append(item)
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
    headers = ["№/Type",'Kolmogorov True', 'Kolmogorov False', 'NIST True', 'NIST False']
    all_seq = []
    gb = []
    with open(name_report_file, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
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
            all_seq.append([t_seq_kolmogorov_count, f_seq_kolmogorov_count, t_seq_nist_count, f_seq_nist_count])
            ii+=1
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
        gst = ['1-50']
        bst = ['51-100']
        for g, b in zip(gb_counter(good_seq_tests), gb_counter(bad_seq_tests)):
            gst.append(g)
            bst.append(b)
        gb.append(gst)
        gb.append(bst)
        writer.writerow(gst)
        writer.writerow(bst)
        return all_seq, gb

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

def task(stat, a, β, length):
    headers = [["alpha", "beta", "bits length"], ["№/Type",'Kolmogorov', 'p_n', 'p_v', "NIST", 'p_n', 'p_v']]
    with open("log.csv", "a", encoding='UTF8', newline='') as logcsv:
        writer = csv.writer(logcsv)
        
        with open("log.txt", "a") as log:   
            if len(stat) != 2:
                ii = 1
                writer.writerow(headers[0])
                writer.writerow([a, β, length])
                writer.writerow('\n')
                writer.writerow(headers[1])
                for s in stat:
                    p1_kol, p2_kol = cipbmie(6, s[1], β)
                    p1_nist, p2_nist = cipbmie(6, s[3], β)
                    writer.writerow([ii, p1_kol <= a <= p2_kol, p1_kol, p2_kol, p1_nist <= a <= p2_nist, p1_nist, p2_nist])
                    
                    if p1_kol <= a <= p2_kol:
                        log.write(f"Kolmogorov: The sequence {ii} is accepted with the significance level {a} and the confidence level {β}. [{p1_kol}; {p2_kol}]\n")
                    else:
                        log.write(f"Kolmogorov: The sequence {ii} is rejected. [{p1_kol}; {p2_kol}]\n")
                    if p1_nist <= a <= p2_nist:
                        log.write(f"NIST: The sequence {ii} is accepted with the significance level {a} and the confidence level {β}. [{p1_nist}; {p2_nist}]\n")
                    else:
                        log.write(f"NIST: The sequence {ii} is rejected. [{p1_nist}; {p2_nist}]\n")
                    ii+=1
            else:
                p1_good_kol, p2_good_kol = cipbmie(300, stat[0][1], β)
                p1_good_nist, p2_good_nist = cipbmie(300, stat[0][3], β)

                p1_bad_kol, p2_bad_kol = cipbmie(300, stat[1][1], β)
                p1_bad_nist, p2_bad_nist = cipbmie(300, stat[1][3], β)
                writer.writerow(["GOOD", p1_good_kol <= a <= p2_good_kol, p1_good_kol, p2_good_kol, p1_good_nist <= a <= p2_good_nist, p1_good_nist, p2_good_nist])
                writer.writerow(["BAD", p1_bad_kol <= a <= p2_bad_kol, p1_bad_kol, p2_bad_kol, p1_bad_nist <= a <= p2_bad_nist, p1_bad_nist, p2_bad_nist])
                writer.writerow('\n')
                if p1_good_kol <= a <= p2_good_kol:
                    log.write(f"GOOD:\n\tKolmogorov: accepted with the significance level {a} and the confidence level {β}. [{p1_good_kol}; {p2_good_kol}]\n")
                else:
                    log.write(f"GOOD:\n\tKolmogorov: rejected. [{p1_good_kol}; {p2_good_kol}]\n")

                if p1_good_nist <= a <= p2_good_nist:
                    log.write(f"GOOD:\n\tNIST: accepted with the significance level {a} and the confidence level {β}. [{p1_good_nist}; {p2_good_nist}]\n")
                else:
                    log.write(f"GOOD:\n\tNIST: rejected. [{p1_good_nist}; {p2_good_nist}]\n")

                if p1_bad_kol <= a <= p2_bad_kol:
                    log.write(f"BAD:\n\tKolmogorov: accepted with the significance level {a} and the confidence level {β}. [{p1_bad_kol}; {p2_bad_kol}]\n")
                else:
                    log.write(f"BAD:\n\tKolmogorov: rejected. [{p1_bad_kol}; {p2_bad_kol}]\n")

                if p1_bad_nist <= a <= p2_bad_nist:
                    log.write(f"BAD:\n\tNIST: accepted with the significance level {a} and the confidence level {β}. [{p1_bad_nist}; {p2_bad_nist}]\n")
                else:
                    log.write(f"BAD:\n\tNIST: rejected. [{p1_bad_nist}; {p2_bad_nist}]\n")

def main(bits_len = [256, 512], num_of_seq = 100, alpha = [0.01, 0.05, 0.1], beta = [0.9, 0.95]):
        clear_folder_from_file_type('xlsx')
        clear_folder_from_file_type('txt')
        log = open("log.txt", 'w')
        log.close()
        for bl in bits_len:
            log = open("log.txt", 'a')
            log.write(f'Length of {bl} bits\n')
            log.close()
            seq = generate_seqs(bl, num_of_seq)
            for a in alpha:
                log = open("log.txt", 'a')
                log.write(f'Significance level {a}\n')
                log.close()
                resume = methods_fill(seq, a)
                res = make_csv(resume, f'Результати тестів послідовності довжиною {bl} біт({int(a*100)}%).csv')
                all_seq_stat, gb_stats = report(res, f'{bl}.{int(a*100)}.stat.csv')
                for b in beta:
                    log = open("log.txt", 'a')
                    log.write(f'Beta {b}\n')
                    log.close()
                    task(all_seq_stat, a, b, bl)
                    task(gb_stats, a, b, bl)
        convert_csv_to_xlsx()
        time.sleep(3)
        clear_folder_from_file_type('csv')
  
if __name__ == "__main__": 
    main()