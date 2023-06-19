import csv

import pandas as pd
from plotnine import *
from utils import *


# get a csv file that contains the shadow region information, the input are annotation files of mature miRNAs
def get_shadow_interval(file):
    with open(file, 'r') as g_f:
        with open('miRNA_data.csv', 'w') as w:
            writer = csv.writer(w)
            for line in g_f:
                if '#' in line:
                    continue
                g_line_list = line.split('\t')
                if g_line_list[2] == 'miRNA_primary_transcript':
                    # p = re.compile('Name=(ath-MIR[0-9]*[a-z]?)')
                    p = re.compile('Name=(.*);?')
                    ID = p.findall(g_line_list[8])[0]
                    if g_line_list[6] == '+':
                        start = int(g_line_list[3].strip())
                    else:
                        start = int(g_line_list[4].strip())

                if g_line_list[2] == 'miRNA':
                    if g_line_list[6] == '+':
                        start1 = int(g_line_list[3].strip()) - start + 1
                        end1 = int(g_line_list[4].strip()) - start + 1

                    else:
                        start1 = start - int(g_line_list[4].strip()) + 1
                        end1 = start - int(g_line_list[3].strip()) + 1
                    l1 = [ID, start1, end1]
                    writer.writerow(l1)


# get a csv file that contains the reads number of, input the gff3 file
def get_reads_number(file, bed_dic):
    with open(file, 'r') as g_f:
        with open('miRNA_reads.csv', 'w') as mi:

            writer = csv.writer(mi)

            for line in g_f:
                if '#' in line:
                    continue
                g_line_list = line.split('\t')

                if g_line_list[2] == 'miRNA_primary_transcript':
                    gff_dic = {}
                    l1 = ['\t']
                    l2 = []
                    # p = re.compile('Name=(ath-MIR[0-9]*[a-z]?)')
                    p = re.compile('Name=(.*);?')
                    p1 = re.compile(r"chr", re.I)
                    s = p1.search(g_line_list[0]).group()
                    ID = p.findall(g_line_list[8])[0]
                    l2.append(ID)

                    if g_line_list[6] == '+':
                        start = int(g_line_list[3].strip())
                        end = int(g_line_list[4].strip())

                        for tmp in range(start, end + 1):
                            k = g_line_list[0].strip(s) + '_' + str(tmp) + '_' + g_line_list[6]
                            gff_dic[str(tmp)] = bed_dic.get(k, 0)
                            l1.append(str(tmp - start + 1))
                            l2.append(gff_dic[str(tmp)])

                    else:
                        start = int(g_line_list[4].strip())
                        end = int(g_line_list[3].strip())

                        for tmp in range(-start, -end + 1):
                            k = g_line_list[0].strip(s) + '_' + str(-tmp) + '_' + g_line_list[6]
                            gff_dic[str(-tmp)] = bed_dic.get(k, 0)
                            l1.append(str(tmp + start + 1))
                            l2.append(gff_dic[str(-tmp)])
                    writer.writerow(l1)
                    writer.writerow(l2)


# get the plot that contains the mirna sites and a shaded section showing the miRNA region
def line_plotting(name_of_mirna):
    X = []
    Y = []
    with open('miRNA_reads.csv', 'r') as mi:
        for row in mi:
            if '\t' in row:
                index_list = row.split(',')
            else:
                if name_of_mirna in row:
                    X = index_list[1::]
                    Y = row.split(',')[1::]
                    break
    if len(X) == 0:
        print('There is a wrong in the name of miRNA!')
        exit()

    start2 = end2 = start1 = end1 = 0
    with open('miRNA_data.csv', 'r') as da:
        repeat = 1
        for line in da:
            line = line.strip()
            if name_of_mirna in line:
                da_list = line.split(',')
                if repeat == 3:
                    break
                if repeat == 2:
                    start2 = int(da_list[1])
                    end2 = int(da_list[2])
                    repeat += 1
                if repeat == 1:
                    start1 = int(da_list[1])
                    end1 = int(da_list[2])
                    repeat += 1

    X[-1] = X[-1].strip()
    Y[-1] = Y[-1].strip()
    label = []
    for i in range(0, len(Y)):
        Y[i] = int(Y[i])
        label.append(name_of_mirna)
    data_info_long = pd.DataFrame({'x': X, 'y': Y, 'label': label})
    data_info_long['x'] = pd.Categorical(data_info_long.x, categories=pd.unique(data_info_long.x))

    region1 = (start1, end1)
    region2 = (start2, end2)
    p = (
            ggplot(data_info_long, aes('x')) +
            geom_line(aes(x='x', y='y', group='label', color='label')) +
            geom_density() +
            annotate(geom_rect, xmin=region1[0], xmax=region1[1], ymin=0, ymax=float('inf'), alpha=0.5) +
            annotate(geom_vline, xintercept=region1, alpha=0) +
            annotate(geom_rect, xmin=region2[0], xmax=region2[1], ymin=0, ymax=float('inf'), alpha=0.5) +  # new line
            annotate(geom_vline, xintercept=region2, alpha=0) +
            labs(x='relative position of bp', y='the number of reads') +
            scale_x_discrete(breaks=[str(x) for x in range(0, len(X), 5)]) +
            theme(axis_text_x=element_text(size=8, family="Monospace", color="black", angle=90),
                  axis_title=element_text(size=12))

    )
    p.save(name_of_mirna + "_analysis_result.pdf")


# print the miRNAs that have an enrichment at 5p or 3p+1
def enrichment_signals(bed_dic, gff_file):
    with open('cleavage_site.csv', 'w') as cl:
        with open(gff_file, 'r') as g_f:
            writer = csv.writer(cl)
            writer.writerow(['Name', '5p', '3p+1', '*5p', '*3p+1'])
            l = []
            for line in g_f:
                if '#' in line:
                    continue
                g_line_list = line.split('\t')
                p1 = re.compile(r"chr", re.I)
                s = p1.search(g_line_list[0]).group()

                if g_line_list[2] == 'miRNA_primary_transcript':
                    # p = re.compile('Name=(ath-MIR[0-9]*[a-z]?)')
                    p = re.compile('Name=(.*);?')
                    ID = p.findall(g_line_list[8])[0]
                    if l[1:]:
                        if max(l[1:]) != 0:
                            writer.writerow(l)
                    l = [ID]
                    Sum = 0
                    length = 0

                    if g_line_list[6] == '+':
                        start = int(g_line_list[3].strip())
                        end = int(g_line_list[4].strip())

                        for tmp in range(start, end + 1):
                            k = g_line_list[0].strip(s) + '_' + str(tmp) + '_' + g_line_list[6]
                            Sum += bed_dic.get(k, 0)
                            length += 1

                    else:
                        start = int(g_line_list[4].strip())
                        end = int(g_line_list[3].strip())

                        for tmp in range(-start, -end + 1):
                            k = g_line_list[0].strip(s) + '_' + str(-tmp) + '_' + g_line_list[6]
                            Sum += bed_dic.get(k, 0)
                            length += 1
                if g_line_list[2] == 'miRNA':

                    if g_line_list[6] == '+':
                        start = g_line_list[3].strip()
                        end = int(g_line_list[4].strip()) + 1
                    else:
                        start = g_line_list[4].strip()
                        end = int(g_line_list[3].strip()) - 1
                    k1 = g_line_list[0].strip(s) + '_' + start + '_' + g_line_list[6]
                    k2 = g_line_list[0].strip(s) + '_' + str(end) + '_' + g_line_list[6]
                    if Sum == 0:
                        continue
                    s_5p = bed_dic.get(k1, 0) / (Sum / length)
                    s_3p = bed_dic.get(k2, 0) / (Sum / length)
                    l.append(s_5p)
                    l.append(s_3p)


def mode3_analysis(gff3_file, name_of_mirna, bed_list):
    bed_dic = get_bed_dicts(bed_list)
    get_shadow_interval(gff3_file)
    get_reads_number(gff3_file, bed_dic)
    line_plotting(name_of_mirna)
    enrichment_signals(bed_dic, gff3_file)
