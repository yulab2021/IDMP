import re
import os


# get the reads information from the bed file(Separate samples)
def get_bed_dict(file):
    bed_dic = {}
    file = file.strip()
    with open(file, 'r') as b_f:
        for line in b_f:
            line = line.strip()
            b_line_list = line.split('\t')
            #flag = b_line_list[6]
            """if flag not in ["0", "16"]:
                continue"""
            if re.search('[1-9]', b_line_list[0]):
                if b_line_list[5] == '+':
                    pos = int(b_line_list[1].strip()) + 1

                else:
                    pos = int(b_line_list[2].strip())

                b_k = b_line_list[0].strip("Chr") + '_' + str(pos) + '_' + b_line_list[5]
                if b_k not in bed_dic:
                    bed_dic[b_k] = 0
                bed_dic[b_k] += 1
    return bed_dic


# Read the information of *.bed file and create a dict
def get_bed_dicts(file_list):
    bed_dic = {}
    with open(file_list, 'r') as f:
        for file in f:
            file = file.strip()
            with open(file, 'r') as b_f:
                for line in b_f:
                    line = line.strip()
                    b_line_list = line.split('\t')

                    # flag = b_line_list[6]
                    """if flag not in ["0", "16"]:
                        continue"""
                    if re.search('[1-9]', b_line_list[0]):
                        if b_line_list[5] == '+':
                            pos = int(b_line_list[1].strip()) + 1

                        else:
                            pos = int(b_line_list[2].strip())

                        b_k = b_line_list[0].strip('Chr') + '_' + str(pos) + '_' + b_line_list[5]
                        if b_k not in bed_dic:
                            bed_dic[b_k] = 0
                        bed_dic[b_k] += 1
    return bed_dic


# Extract the reads information of the CDS region using the comment file
def get_gff_dict_cds(file, bed_dic):
    with open(file, 'r') as g_f:

        gff_dic = {}
        gff_dic1 = {}
        id_list = []

        for line in g_f:
            if '#' in line:
                continue
            line = line.strip()
            g_line_list = line.split('\t')
            if re.search('[0-9]', g_line_list[0]):
                if g_line_list[2] == 'CDS':
                    if re.search('Parent=(.*\.1),?', g_line_list[-1]):
                        p = re.compile('Parent=([^,]*\.1),?')
                        ID = p.findall(g_line_list[8])[0]
                        if "," in ID:
                            ID = ID.split(",")[0]
                        if ID not in id_list:
                            id_list.append(ID)
                            flag = 1
                        else:
                            flag = 0

                        if flag:
                            i = 1
                            gff_dic1[ID] = 0
                        strand = g_line_list[6]
                        if strand == '+':
                            start = int(g_line_list[3].strip())
                            end = int(g_line_list[4].strip())
                            gff_dic1[ID] += end - start + 1
                            for tmp in range(start, end + 1):
                                k = g_line_list[0].strip("Chr") + '_' + str(tmp) + '_' + strand
                                gff_dic[ID + '_' + str(i)] = bed_dic.get(k, 0)
                                i += 1

                        else:
                            start = int(g_line_list[4].strip())
                            end = int(g_line_list[3].strip())
                            gff_dic1[ID] += start - end + 1
                            for tmp in range(-start, -end + 1):
                                k = g_line_list[0].strip("Chr") + '_' + str(-tmp) + '_' + strand
                                gff_dic[ID + '_' + str(i)] = bed_dic.get(k, 0)
                                i += 1

    return gff_dic, gff_dic1


# prepare a dict to record the codons and their reads number
def get_codon_reads_dict():
    codon_list = []
    codon_reads_dic = {}
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            for ch in ['A', 'C', 'G', 'T']:
                s = i + j + ch
                codon_list.append(s)
                codon_reads_dic[s] = 0
    return codon_list, codon_reads_dic


# Obtain a dictionary containing the sequence of bases corresponding to the segment
def get_reference_dict(file):
    file = file.strip()
    with open(file) as f:
        ge_dict = {}
        sequence = ''
        for line in f:
            line = line.strip()
            if '>' in line:
                p = re.compile('>[Cc]?h?r?([0-9]+)')
                if len(p.findall(line)) == 0:
                    ge_dict[ID] = sequence
                    break
                else:
                    if sequence != '':
                        ge_dict[ID] = sequence
                    sequence = ''
                    ID = str(p.findall(line)[0])
            else:
                sequence += line
    return ge_dict


# get the sequence of cds, the file 'cds.fa' can be used to Submit to the site
def extract_cds(genome_file, gff_file):
    cmd = f"gffread {gff_file} -g {genome_file} -x cds.fa"
    os.system(cmd)


# get a dict that contains the gene id, and it's corresponding gene name
def get_gene_name_dict(gene_name_file):
    gene_name_dic = {}
    with open(gene_name_file, 'r') as f:
        for line in f:
            line_list = line.split('\t')
            try:
                if line_list[-2].strip() != '':
                    gene_name_dic[line_list[1].strip() + '.1'] = line_list[-2].strip()
                else:
                    gene_name_dic[line_list[1].strip() + '.1'] = line_list[1].strip() + '.1'

            except IndexError:
                continue
        """p = re.compile('AT[0-5]G[0-9]+')
        for line in f:
            if p.search(line) is None:
                continue
            line_list = line.strip().split('\t')
            gene_name_dic[line_list[0].strip()] = line_list[1].strip()"""
    return gene_name_dic

# 'Parent=(AT[0-5]G[0-9]+\.1)'


def mysum(l):
    Sum = 0
    for i in l:
        Sum += int(i)
    return Sum

