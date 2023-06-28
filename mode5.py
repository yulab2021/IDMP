import csv
import pandas as pd
import seaborn as sns
import numpy as np
from utils import *
import matplotlib.pyplot as plt
from plotnine import *
from sklearn import preprocessing


# get a csv file that contains EJC reads information and TSI
def get_ejc_csv_file(file_name, bed_dic, gff_file, gene_name_file):
    gene_name_dict = get_gene_name_dict(gene_name_file)
    with open(gff_file, 'r') as g_f:
        with open(file_name.strip("Aligned.sortedByCoord.out.bed") + "_EJC_upstream_50.csv", 'w') as w:

            writer = csv.writer(w)
            l1 = list()
            l1.append('gene id')
            for i in range(-50, 1):
                l1.append(str(i))
            l1.append('TSI')
            l1.append("gene name")
            writer.writerow(l1)

            id_list = []
            num_of_repeat = 0

            for line in g_f:
                if '#' in line:
                    continue
                line = line.strip()
                g_line_list = line.split('\t')
                if g_line_list[2] == 'exon':
                    # if re.search('Parent=(AT[0-9]G[0-9])', g_line_list[8]):
                    if re.search('Parent=(.*),?', g_line_list[8]):
                        p = re.compile('Parent=([^,]*\.1),?')
                        ID_list = p.findall(g_line_list[8])
                        if ID_list:
                            ID = ID_list[0]
                        else:
                            continue
                        if ID not in id_list:
                            id_list.append(ID)
                            num_of_repeat = 1
                            ID += '_' + str(1)
                        else:
                            ID += '_' + str(num_of_repeat + 1)
                            num_of_repeat += 1
                        l2 = list()
                        l2.append(ID)
                        total = 0
                        if g_line_list[6] == '+':
                            ejc = int(g_line_list[4].strip())
                            for tmp in range(ejc - 50, ejc + 1):
                                k = g_line_list[0].strip('Chr') + '_' + str(tmp) + '_' + g_line_list[6]
                                l2.append(bed_dic.get(k, 0))
                                total += bed_dic.get(k, 0)
                        else:
                            ejc = int(g_line_list[3].strip())
                            ejc = -ejc
                            for tmp in range(ejc - 50, ejc + 1):
                                k = g_line_list[0].strip('Chr') + '_' + str(-tmp) + '_' + g_line_list[6]
                                l2.append(bed_dic.get(k, 0))
                                total += bed_dic.get(k, 0)

                        ave = int(l2[23]) + int(l2[24]) + int(l2[25])
                        ave /= 3
                        # ave /= (mysum(l2[1:]) / len(l2[1:]))
                        if total == 0:
                            ave = 0
                        else:
                            ave /= total / len(l2[1:])
                        if 2 < ave < 50:
                            l2.append(ave)
                            if max(l2[1:23]) < 1000 and max(l2[26:-1]) < 1000:
                                l2.append(gene_name_dict.get(ID, ID))
                                writer.writerow(l2)


# the file_list contains all the bed files, the sra_file is an
# annotate document contains the temp and tissue information of the samples
def get_file_for_plotting(file_list, gff_file, sra_file):
    with open("EJC_upstream_50.csv", 'w') as w:
        writer = csv.writer(w)
        l1 = list()
        l1.append('\t')
        for i in range(-50, 1):
            l1.append(str(i))
        l1.append('Treatment')
        l1.append('Condition')
        writer.writerow(l1)

    xlsx_list = []
    with open(sra_file, 'r') as f:
        for line in f:
            if "Run" in line:
                continue
            temp = line.split(',')[-2]
            tissue = line.split(',')[-1]
            xlsx_list.append(temp)
            xlsx_list.append(tissue)

    num_of_row = 0  # Used when adding temperature and tissue to the list
    with open(file_list, 'r') as fl:
        with open("EJC_upstream_50.csv", 'a') as w:
            writer = csv.writer(w)
            for file_name in fl:
                file_name = file_name.strip()
                l2 = [0] * 52
                l2[0] = file_name.strip('Aligned.sortedByCoord.out.bed')
                bed_dic = get_bed_dict(file_name)
                """with open(file_name, 'r') as b_f:
                    bed_dic = {}
                    for line in b_f:
                        line = line.strip()
                        b_line_list = line.split('\t')
                        if re.search('[1-9]', b_line_list[0]):
                            if b_line_list[5] == '+':
                                pos = int(b_line_list[1].strip()) + 1
                                x = b_line_list[0] + '_' + str(pos) + '_' + b_line_list[5]

                            else:
                                pos = int(b_line_list[2].strip())
                                x = b_line_list[0] + '_' + str(pos) + '_' + b_line_list[5]
                            if x not in bed_dic.keys():
                                bed_dic[x] = 0
                            bed_dic[x] += 1"""

                with open(gff_file, 'r') as g_f:

                    for line in g_f:
                        if '#' in line:
                            continue
                        line = line.strip()
                        g_line_list = line.split('\t')
                        if g_line_list[2] == 'exon':
                            # if re.search('Parent=AT[0-5]G[0-9]+\.1', g_line_list[8]):
                            if re.search('Parent=(.+\.1),?', g_line_list[8]):
                                if g_line_list[6] == '+':
                                    ejc = int(g_line_list[4].strip())

                                    for tmp in range(ejc - 50, ejc + 1):
                                        k = g_line_list[0].strip('Chr') + '_' + str(tmp) + '_' + g_line_list[6]
                                        l2[tmp + 50 - ejc + 1] += bed_dic.get(k, 0)

                                else:
                                    ejc = int(g_line_list[3].strip())
                                    ejc = -ejc

                                    for tmp in range(ejc - 50, ejc + 1):
                                        k = g_line_list[0].strip('Chr') + '_' + str(-tmp) + '_' + g_line_list[6]

                                        l2[tmp + 50 - ejc + 1] += bed_dic.get(k, 0)

                sum_row = sum(l2[1:52:])
                for i in range(51):
                    l2[i + 1] = l2[i + 1] / sum_row
                    l2[i + 1] *= 100

                l2.append(xlsx_list[num_of_row])
                l2.append(xlsx_list[num_of_row + 1].strip('\n'))
                num_of_row += 2

                writer.writerow(l2)


# file_name: a csv file used to plot
def ejc_plotting(file_name):
    data_info = pd.read_csv(file_name)
    legend_list = []
    repeat_list = []
    with open(file_name, 'r') as df:
        for line in df:
            if str(-50) in line:
                continue
            line_list = line.split(',')
            tag = line_list[-2] + '_' + line_list[-1].strip("\"").strip('\n')
            if tag not in repeat_list:
                repeat_list.append(tag)
                repeat = 1
            else:
                repeat += 1

            legend_list.append(tag + '_' + str(repeat))

    data_info.insert(data_info.shape[1], 'Sample', legend_list)

    data_info = data_info.T

    # wide table to long table

    X = []
    Y = []
    sample = []

    index_list = list(data_info.index)
    k = 1
    for i in range(1, 51):
        loc = np.array(data_info.iloc[i])
        for j in range(len(loc)):
            Y.append(loc[j])
            X.append(index_list[k])
            sample.append(legend_list[j])
        k += 1

    data_info_long = pd.DataFrame({"x": X, "y": Y, "sample": sample})
    data_info_long['x'] = pd.Categorical(data_info_long.x, categories=pd.unique(data_info_long.x))

    p = (
            ggplot(data_info_long) +
            geom_line(aes(x='x', y='y', group="sample", color="sample"), show_legend=True) +  # group
            labs(x='Distance from exon 3\' end (nt)', y='Relative frequency of 5\'P end occurrences (%)') +

            theme(legend_position=(0.25, 0.9), legend_text=element_text(size=6), legend_key_size=12,
                  legend_background=element_rect(alpha=0), figure_size=(12, 8), axis_text=element_text(size=8),
                  axis_title=element_text(size=12))

    )

    p.save(file_name.strip('.csv') + '_' + ".pdf")


# input the csv file that contains TSI
def get_ejc_heatmap(file_name):
    data = pd.read_csv(file_name.strip("Aligned.sortedByCoord.out.bed") + "_EJC_upstream_50.csv", index_col=0)
    data = data.sort_values(by='TSI', ascending=False)
    data.pop("TSI")
    data.pop("gene name")
    scaler = preprocessing.StandardScaler().fit(data.T)  # Rows were normalized
    data_T_scale = scaler.transform(data.T)
    data_scale = data_T_scale.transpose()
    plt.figure()
    row_labels = list(data.index)
    col_labels = list(data.columns)
    df = pd.DataFrame(data_scale.T, index=col_labels, columns=row_labels)
    sns.clustermap(df.T, row_cluster=False, col_cluster=False, cmap="OrRd", center=0, yticklabels='')
    # plt.xlim(-50, 0)
    # Spectral coolwarm
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(file_name.strip("Aligned.sortedByCoord.out.bed") + "_ejc_heatmap.pdf")
