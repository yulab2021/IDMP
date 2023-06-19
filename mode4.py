import csv
import math
import pandas as pd
import numpy as np
from plotnine import *
from utils import *


# get the sequence of miRNA, the file 'mirna.fa' can be used to Submit to the site
# the mature_file is the mature RNA annotation files, which a can be downloaded from www.mirbase.org
def extract_mirna(mature_file, rna_id):
    with open(mature_file, 'r') as fa:
        with open("mirna.fa", 'w') as mi:
            flag = False
            for line in fa:
                if flag:
                    mi.writelines(line)
                    flag = False
                r_id = rna_id + '-miR'
                if r_id in line:
                    mi.writelines(line)
                    flag = True


# get a special cds_dict that contains the reads of upstream and downstream 50 nt
def get_cds_dict(file, bed_dic):
    cds_dict = {}
    with open(file, 'r') as g_f:
        id_list = []
        strand = ''
        for row in g_f:
            if "#" in row:
                continue
            row_list = row.split('\t')
            if row_list[2] == 'CDS':
                # if re.search('Parent=(AT[0-5]G[0-9]+\.1)', row_list[8]):
                # p = re.compile('Parent=(AT[0-5]G[0-9]+\.1)')
                if re.search('Parent=(.*\.1),?', row_list[-1]):
                    p = re.compile('Parent=(.*\.1),?')
                    ID = p.findall(row_list[-1])[0]
                    if ID not in id_list:
                        try:
                            pre_id = id_list[-1]
                        except IndexError:
                            pre_id = ID
                        if strand:
                            if strand == "+":
                                for tmp in range(end + 1, end + 51):
                                    k = Chr + '_' + str(tmp) + '_' + strand
                                    cds_dict[pre_id + '_' + str(i)] = bed_dic.get(k, 0)
                                    i += 1
                            else:
                                for tmp in range(-end + 1, -end + 51):
                                    k = Chr + '_' + str(-tmp) + '_' + strand
                                    cds_dict[pre_id + '_' + str(i)] = bed_dic.get(k, 0)
                                    i += 1

                        id_list.append(ID)
                        flag = 1
                    else:
                        flag = 0
                        d = 0

                    if flag:
                        i = -50
                        d = 50

                    strand = row_list[6]
                    Chr = row_list[0].strip('Chr')
                    if strand == "+":
                        start = int(row_list[3])
                        end = int(row_list[4])

                        for tmp in range(start - d, end + 1):
                            k = Chr + '_' + str(tmp) + '_' + strand
                            cds_dict[ID + '_' + str(i)] = bed_dic.get(k, 0)
                            i += 1

                    else:
                        start = int(row_list[4])
                        end = int(row_list[3])

                        for tmp in range(-start - 50, -end + 1):
                            k = Chr + '_' + str(-tmp) + '_' + strand
                            cds_dict[ID + '_' + str(i)] = bed_dic.get(k, 0)
                            i += 1

    return cds_dict


# the result_file is the file downloaded from the website. We process the file and add some new information
def mode4_analysis(file_list, gff_file_name, result_file):
    bed_dic = get_bed_dicts(file_list)
    cds_dic = get_cds_dict(gff_file_name, bed_dic)
    with open(result_file, 'r') as t_file:
        with open("mirna_5p_reads.csv", 'w') as out:
            writer = csv.writer(out)
            id_list = []
            short_id_list = []
            for line in t_file:
                if '#' in line:
                    continue
                if "miRNA" in line:
                    title_list = line.split('\t')
                    title_list.extend(["miRNA cleavage site signal", "mean flanking signal"
                                          , "miRNA cleavage efficiency"])
                    writer.writerow(title_list)

                # p = re.compile('AT[0-5]G[0-9]+\.1')
                # if p.findall(line):

                line_list = line.split('\t')
                if line_list[11] == "Cleavage":

                    start = int(line_list[6])
                    end = int(line_list[7])
                    ID = line_list[1].split('.')[0] + '.' + line_list[1].split('.')[1]
                    p = re.compile('\.1{1}$')
                    if not p.search(ID):
                        continue
                    p1 = re.compile('.*-miR[1-9]+')
                    long_id = line_list[0]
                    short_id = p1.findall(line_list[0])[0]

                    if long_id not in id_list and short_id in short_id_list:
                        continue
                    if long_id not in id_list and short_id not in short_id_list:
                        id_list.append(long_id)
                        short_id_list.append(short_id)

                    Site = end - 10
                    # Sum = 0
                    data_l = []
                    for tmp in range(Site - 50, Site + 51):
                        if tmp == Site:
                            continue
                        else:
                            data_l.append(cds_dic[ID + '_' + str(tmp)])
                            # Sum += cds_dic[ID + '_' + str(tmp)]
                    site_reads = cds_dic[ID + '_' + str(Site)]
                    if site_reads < 5:
                        continue
                    medium = np.median([i for i in data_l if i != 0])
                    # Sum = Sum / float(100)
                    try:
                        efficiency = math.log2(site_reads / medium)
                    except Exception:
                        efficiency = 0
                    if efficiency >= 1:
                        line_list[-1] = line_list[-1].strip()
                        line_list.extend([site_reads, medium, efficiency])
                        writer.writerow(line_list)


def reads_plotting(mirna_name, bed_files, gff_file):
    site_dic = {}
    bed_dic = get_bed_dicts(bed_files)
    gff_dic, gff_dic1 = get_gff_dict_cds(gff_file, bed_dic)
    with open("mirna_5p_reads.csv", 'r') as f:
        for line in f:
            line_list = line.split(',')
            if line_list[0].strip("\"").strip() == mirna_name:
                ID = line_list[1].split(".")[0] + ".1"
                for i in range(-50, 51):
                    site_dic[i] = gff_dic.get((ID + "_" + str(int(line_list[6]) + i)), 0)
                plotting(site_dic, ID)
                site_dic = {}


def plotting(data: dict, name_of_gene):
    col1 = []
    col2 = []
    for k in data:
        col1.append(k)
        col2.append(data[k])
    data_info_long = pd.DataFrame({"x": col1, "y": col2, "label": name_of_gene})
    data_info_long['x'] = pd.Categorical(data_info_long.x, categories=pd.unique(data_info_long.x))

    p = (
            ggplot(data_info_long) +
            geom_line(aes(x='x', y='y', group='label', color='label')) +
            labs(x='Number relative to miRNA start site(nt)', y='Read count') +
            scale_x_discrete(breaks=[5 * i for i in range(-10, 11)]) +
            theme(legend_position=(0.25, 0.9), legend_box='vertical', legend_text=element_text(size=6),
                  legend_key_size=12, legend_margin=1.5, legend_background=element_rect(alpha=0),
                  figure_size=(12, 8), axis_text=element_text(size=8), axis_title=element_text(size=12))
    )
    p.save(name_of_gene + "_reads_flanking_100.pdf")
