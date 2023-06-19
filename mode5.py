import csv
import pandas as pd
import seaborn as sns
import numpy as np
from sklearn import preprocessing
from utils import *
import matplotlib.pyplot as plt


def uorf_analysis(file_names, gff_file, gene_name_dic, genome_file):
    # get five_prime_UTR sequence
    with open(file_names, 'r') as f:
        for file in f:
            file = file.strip()
            with open(gff_file, 'r') as g_f:

                id_list = []
                bed_dic = get_bed_dict(file)
                ge_dic = get_reference_dict(genome_file)
                utr_dict = {}
                i = 1
                for line in g_f:
                    if "#" in line:
                        continue
                    line_list = line.split('\t')
                    p1 = re.compile('[0-9]+')
                    p = re.compile('Parent=(.*\.1)\.?')
                    if len(p1.findall(line_list[0])) != 0 and len(p.findall(line_list[-1])) != 0:
                        if line_list[2] == "five_prime_UTR":

                            Chr = line_list[0].strip('Chr')
                            # ID = line.split("=")[1].strip()
                            ID = p.findall(line_list[-1])[0]

                            if ID in id_list:
                                ID += f'-{i}'
                                i += 1
                            else:
                                i = 1
                            id_list.append(ID)
                            start = int(line_list[3])
                            end = int(line_list[4])
                            strand = line_list[6]
                            utr_dict[Chr + '_' + ID + '_' + strand] = [start, end]
            se_dict = {}
            reads_dict = {}
            for key in utr_dict:
                strand = key.split('_')[-1]
                p2 = re.compile('\d_(.*)_+|-')
                ID = p2.findall(key)[0]
                Chr = key.split('_')[0]
                start = utr_dict[key][0]
                end = utr_dict[key][1]
                if strand == '+':
                    se_dict[ID] = ge_dic[Chr][start:end + 1]
                    for i in range(start - 50, end + 51):
                        k0 = "%s_%s" % (ID, str(i - start))
                        k = "%s_%s_%s" % (Chr, i + 1, strand)
                        reads_dict[k0] = bed_dic.get(k, 0)
                else:
                    se_dict[ID] = ge_dic[Chr][::-1][start:end + 1]
                    for i in range(start - 50, end + 51):
                        k0 = "%s_%s" % (ID, str(i - start))
                        k = "%s_%s_%s" % (Chr, len(ge_dic[Chr]) - i, strand)
                        reads_dict[k0] = bed_dic.get(k, 0)
            # get a file that contains the upstream reads number of 5'uorf
            with open(file.strip("Aligned.sortedByCoord.out.bed") + 'uorf.csv', 'w') as out:
                writer = csv.writer(out)
                title = ['\t', 'gene name']
                for i in range(-48, 3):
                    title.append(i)
                title.append("TSI")
                writer.writerow(title)
                for ID in id_list:
                    length = len(se_dict[ID])
                    c_list = []
                    for j in range(length - 2):
                        codon = se_dict[ID][j:j + 3]
                        if codon in ["ATG", "CTG"]:
                            c = j
                            while 1:
                                c += 3
                                if c > length:
                                    break
                                codon = se_dict[ID][c:c + 3]
                                if codon in ["TAA", "TGA", "TAG"]:
                                    if c in c_list:
                                        break
                                    else:
                                        c_list.append(c)
                                    """reads_list = [ID + '_' + str(j) + '_' + str(c),
                                                  gene_name_dic.get(ID.split('.')[0], ID)]"""
                                    reads_list = [ID + '_' + str(j) + '_' + str(c),
                                                  gene_name_dic.get(ID, ID)]
                                    for site in range(-48, 3):
                                        k = "%s_%s" % (ID, str(c + site))
                                        reads_list.append(reads_dict[k])
                                    data_l = reads_list.copy()
                                    for site in range(3, 53):
                                        k = "%s_%s" % (ID, str(c + site))
                                        data_l.append(reads_dict.get(k, 0))

                                    """try:
                                        TSI = (int(reads_list[33]) + int(reads_list[34])) / 2 / (
                                                sum(reads_list[2:]) / len(reads_list[2:]))
                                    except ZeroDivisionError:
                                        TSI = 0"""
                                    if not [i for i in reads_list if i != 0]:
                                        if sum(data_l) > 20:
                                            TSI = (int(reads_list[33]) + int(reads_list[34])) / 2 / np.median([i for i in reads_list if i != 0])
                                    else:
                                        TSI = 0
                                    if 3 < TSI < 50 and (int(reads_list[33]) + int(reads_list[34])) >= 2:
                                        reads_list.append(TSI)
                                        writer.writerow(reads_list)
                                    break


def get_uorf_heatmap(file_name):
    file_name = file_name.strip()
    data = pd.read_csv(file_name.strip("Aligned.sortedByCoord.out.bed") + "uorf.csv", index_col=0)
    data = data.sort_values(by='TSI', ascending=False)
    data.pop("TSI")
    data.pop('gene name')
    scaler = preprocessing.StandardScaler().fit(data.T)  # 对行做标准化处理
    data_T_scale = scaler.transform(data.T)
    data_scale = data_T_scale.transpose()
    plt.figure()
    row_labels = list(data.index)
    col_labels = list(data.columns)
    df = pd.DataFrame(data_scale.T, index=col_labels, columns=row_labels)
    sns.clustermap(df.T, row_cluster=False, col_cluster=False, cmap="OrRd", center=0, yticklabels='')
    # Spectral coolwarm
    # plt.xlim(-48, 2)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(file_name.strip("Aligned.sortedByCoord.out.bed") + "uorf_heatmap.pdf")
