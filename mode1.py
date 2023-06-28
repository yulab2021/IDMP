import csv
import numpy as np
from plotnine import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import *
from sklearn import preprocessing
import os


# get the line plot of co_translation
def co_translation_plotting(file_name):
    col1 = []
    col2 = []
    col3 = []

    with open(file_name, 'r', encoding='utf-8') as co:
        reader = csv.reader(co)
        for row in reader:
            col1.append(int(row[0]))
            col2.append(int(row[1]))
            col3.append(row[2])

    data_info_long = pd.DataFrame({"x": col1, "y": col2, "sample": col3})
    data_info_long['x'] = pd.Categorical(data_info_long.x, categories=pd.unique(data_info_long.x))

    p = (
            ggplot(data_info_long) +
            geom_line(aes(x='x', y='y', group="sample", color="sample"), show_legend=True) +  # group
            labs(x='Distance from stop codon(nt)', y='Read count') +
            theme(legend_position=(0.25, 0.9), legend_box='vertical', legend_text=element_text(size=6),
                  legend_key_size=12, legend_margin=1.5, legend_background=element_rect(alpha=0),
                  figure_size=(12, 8), axis_text=element_text(size=8), axis_title=element_text(size=12))

    )

    p.save('c' + file_name.strip('.csv') + ".pdf")


def mode_1_step_1(bed_files, gff_file):
    with open(bed_files, 'r') as f:
        with open('co_translation_plotting.csv', 'w') as co:
            writer = csv.writer(co)
            for file_name in f:
                file_name = file_name.strip()
                bed_dic = get_bed_dict(file_name)
                gff_dic, gff_dic1 = get_gff_dict_cds(gff_file, bed_dic)
                dict1 = {}
                for gene in gff_dic1:
                    end = int(gff_dic1[gene])
                    start = end - 50
                    if start < 0:
                        continue
                    for k in range(start, end + 1):
                        dict1[gene + '_' + str(k - start - 50)] = gff_dic[gene + '_' + str(k)]
                d1_key_list = list(dict1.keys())
                num_list = [0] * 51
                for key in d1_key_list:
                    num = int(key.split('_')[-1])
                    if num != -18 or -19:
                        if int(dict1[key]) > 100:
                            num_list[-num] += 100
                        else:
                            num_list[-num] += int(dict1[key])
                    else:
                        num_list[-num] += int(dict1[key])
                tag = file_name.strip("Aligned.sortedByCoord.out.bed")
                for num in range(-48, 3):
                    w_l = [str(num), str(num_list[-num + 2]), tag]
                    writer.writerow(w_l)

    co_translation_plotting('co_translation_plotting.csv')


def get_ctrd_heatmap(file_name):
    data = pd.read_csv(file_name, index_col=0)
    data = data.sort_values(by='TSI', ascending=False)
    if data.empty:
        return
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
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(file_name.strip('.csv') + "_heatmap.pdf")


# get the csv file of ctrd, which contains the reads of every site and the TSI value
def get_csv_file(file_names, gff_filename, gene_name_file):
    with open(file_names, 'r') as f:
        gene_name_dict = get_gene_name_dict(gene_name_file)
        for file_name in f:
            file_name = file_name.strip('\n')
            bed_dic = get_bed_dict(file_name)
            gff_dic, gff_dic1 = get_gff_dict_cds(gff_filename, bed_dic)
            new_file_name = file_name.strip("Aligned.sortedByCoord.out.bed") + "_co_translation.csv"
            with open(new_file_name, 'w') as w:
                writer = csv.writer(w)
                l2 = ['gene id']
                for i in range(-48, 3):
                    l2.append(str(i))
                l2.append('TSI')
                l2.append("gene name")
                writer.writerow(l2)
                l1 = []
                for gene in gff_dic1:
                    end = int(gff_dic1[gene])
                    start = end - 50
                    if start < 0:
                        continue
                    l1.append(gene)
                    data_l = []
                    for k in range(start, end + 1):
                        l1.append(gff_dic.get(gene + '_' + str(k), 0))
                        data_l.append(gff_dic.get(gene + '_' + str(k), 0))
                    for k in range(end + 1, end + 51):
                        data_l.append(gff_dic[gene + '_' + str(k)])
                    if not [i for i in l1 if i != 0]:
                        if sum(data_l) > 20:
                            TSI = (int(l1[33]) + int(l1[34])) / 2 / np.median([i for i in data_l if i != 0])
                    else:
                        TSI = 0
                    """try:
                        # TSI = (int(l1[33]) + int(l1[34])) / 2 / (mysum(l1[1:]) / len(l1[1:]))
                        TSI = (int(l1[33]) + int(l1[34])) / 2 / np.median([i for i in l1 if i != 0])
                    except ZeroDivisionError:
                        TSI = 0"""
                    if 3 < TSI < 50:
                        l1.append(TSI)
                        l1.append(gene_name_dict[gene])
                        writer.writerow(l1)
                    l1.clear()
            get_ctrd_heatmap(new_file_name)


# Read the information of *.bed file and create a dict for plotting
def pattern_3(file_names, gff_file):
    with open(file_names, 'r') as f:

        new_file_name_list = []
        for file_name in f:
            file_name = file_name.strip('\n')
            bed_dic = get_bed_dict(file_name)
            # Data comparison using the gff file

            gff_dic, gff1_dic = get_gff_dict_cds(gff_file, bed_dic)
            # number of reads over pattern3
            n_gff_dic = sorted(gff_dic)
            new_file_name = file_name.strip("Aligned.sortedByCoord.out.bed") + '_' + "pattern3.csv"
            new_file_name_list.append(new_file_name)
            with open(new_file_name, 'w') as w:
                writer = csv.writer(w)
                l1 = ['gene_id', '0', '1', '2']
                writer.writerow(l1)
                pattern_dict = {}
                l2_dict = {}
                for ID in n_gff_dic:
                    ID_list = ID.split('_')
                    num = int(ID_list[-1]) % 3
                    if ID_list[-2] not in l2_dict:
                        if len(l2_dict) != 0:
                            w_list = [pattern_dict[3], pattern_dict[1], pattern_dict[2], pattern_dict[0]]
                            writer.writerow(w_list)
                        l2_dict[ID_list[-2]] = 1
                        pattern_dict = {3: ID_list[-2]}

                    if ID_list[-1] not in pattern_dict:
                        pattern_dict[num] = int(gff_dic[ID])
                    else:
                        pattern_dict[num] += int(gff_dic[ID])
    return new_file_name_list


# ger a bar plot shows the 3_pattern of reads
def pattern_3_plotting(file_names):
    for file in file_names:
        data_info = pd.read_csv(file)

        X = [0, 1, 2]

        arr0 = pd.to_numeric(data_info['0'])
        arr1 = pd.to_numeric(data_info['1'])
        arr2 = pd.to_numeric(data_info['2'])
        l0 = arr0.tolist()
        l1 = arr1.tolist()
        l2 = arr2.tolist()

        Y = [sum(l0), sum(l1), sum(l2)]
        data_info_long = pd.DataFrame({"x": X, "y": Y})

        p = (
                ggplot(data_info_long) +
                geom_bar(aes(x='x', y='y', fill='x'), stat='identity') +  # group
                labs(x='position in frame', y='reads of 5\'P end') +
                theme(axis_text=element_text(size=8), legend_text=element_text(size=6),
                      axis_title=element_text(size=12))
        )

        p.save(file.strip('.csv') + ".pdf")


# build cds dictionary
def get_cds_dictionary(file):
    with open(file) as f:
        cds_dict = {}  # save genes and corresponding CDS range
        chr_dict = {}  # save genes and corresponding chromosome number
        strand_dict = {}  # save genes and corresponding strand information

        for line in f:
            if "#" in line:
                continue
            type = line.split("\t")[2]
            if type == "CDS":
                chr = line.split('\t')[0].strip('Chr')
                try:
                    int(chr)
                except ValueError:
                    continue
                if re.search('Parent=(.*\.1)', line.split('\t')[-1]):
                    p = re.compile('Parent=([^,]*\.1),?')
                    transcript_id = p.findall(line.split('\t')[-1])[0]
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4])
                strand = line.split("\t")[6]
                if not transcript_id.strip().endswith('.1'):
                    continue
                if transcript_id not in cds_dict:
                    cds_dict[transcript_id] = []
                    chr_dict[transcript_id] = chr
                    strand_dict[transcript_id] = strand
                cds_dict[transcript_id].append([start, end])
    return cds_dict, chr_dict, strand_dict


# from cds.fa extract the sequence odf cds
def get_cds_sequence(file):
    cds_sequence_dict = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if ">" in line:
                transcript_id = line.strip('>').split('.')[0] + '.' + line.split('.')[1]
                cds_sequence_dict[transcript_id] = ""
            else:
                cds_sequence_dict[transcript_id] += line
    return cds_sequence_dict


def transcript_location_to_chr_location(transcript_location, strand, cds_list):
    summation = 0
    for cds in cds_list:
        start = cds[0]
        end = cds[1]
        if transcript_location > summation + end - start + 1:
            summation = summation + end - start + 1
        else:
            if strand == "+":
                chr_location = transcript_location - summation + start - 1
            else:
                chr_location = end - (transcript_location - summation) + 1
            return chr_location


# plot the reads of every codon
def plot_reads_count(codon_reads_dict=None, file_name=None):
    X = []
    Y = []
    for codon in codon_reads_dict:
        if codon_reads_dict[codon] and re.match('[A|C|G|T]{3}', codon):
            X.append(codon)
            Y.append(codon_reads_dict[codon])

    data_dict = dict(zip(X, Y))
    codon_dic1 = dict(sorted(data_dict.items(), key=lambda x: x[1]))
    data_info_long = pd.DataFrame({'codon': codon_dic1.keys(), 'reads': codon_dic1.values()})

    data_info_long['codon'] = pd.Categorical(data_info_long.codon, categories=pd.unique(data_info_long.codon))

    p = (
            ggplot(data_info_long) +
            geom_bar(aes(x='codon', y='reads', fill='codon'), stat='identity') +
            labs(x='codons', y='reads numbers of codon') +
            theme(axis_text_x=element_text(size=8, family="Monospace", color="black", angle=90),
                  figure_size=[8.0, 6.0], axis_title=element_text(size=12))
    )

    p.save(file_name.strip("Aligned.sortedByCoord.out.bed") + '_reads_numbers' + ".pdf")


# plot the counts number of every codon
def plot_counts_count(codon_counts_dict=None, file_name=None):
    X = []
    Y = []
    for codon in codon_counts_dict:
        if codon_counts_dict[codon] and re.match('[A|C|G|T]{3}', codon):
            X.append(codon)
            Y.append(codon_counts_dict[codon])

    data_dict = dict(zip(X, Y))
    codon_dic1 = dict(sorted(data_dict.items(), key=lambda x: x[1]))
    data_info_long = pd.DataFrame({'codon': codon_dic1.keys(), 'counts': codon_dic1.values()})

    data_info_long['codon'] = pd.Categorical(data_info_long.codon, categories=pd.unique(data_info_long.codon))

    p = (
            ggplot(data_info_long) +
            geom_bar(aes(x='codon', y='counts', fill='codon'), stat='identity') +
            labs(x='codons', y='counts numbers of codon') +
            theme(axis_text_x=element_text(size=6.5, family="Monospace", color="black", angle=90),
                  figure_size=[8.0, 6.0], axis_title=element_text(size=12))
    )

    p.save(file_name.strip("Aligned.sortedByCoord.out.bed") + '_reads_numbers' + ".pdf")


# reads number / counts number
def plot_norm_count(codon_reads_dict=None, codon_counts_dict=None, file_name=None):
    X = []
    Y = []
    for codon in codon_reads_dict:
        if codon_reads_dict[codon] and re.match('[A|C|G|T]{3}', codon):
            X.append(codon)
            Y.append(codon_reads_dict[codon] / codon_counts_dict[codon])

    data_dict = dict(zip(X, Y))
    codon_dic1 = dict(sorted(data_dict.items(), key=lambda x: x[1]))
    data_info_long = pd.DataFrame({'codon': codon_dic1.keys(), 'counts': codon_dic1.values()})

    data_info_long['codon'] = pd.Categorical(data_info_long.codon, categories=pd.unique(data_info_long.codon))

    p = (
            ggplot(data_info_long) +
            geom_bar(aes(x='codon', y='counts', fill='codon'), stat='identity') +
            labs(x='codons', y='relative reads of codon') +
            theme(axis_text_x=element_text(size=8, family="Monospace", color="black", angle=90),
                  axis_title=element_text(size=12), figure_size=[8.0, 6.0])
    )

    p.save(file_name.strip("Aligned.sortedByCoord.out.bed") + '_norm_counts' + ".pdf")


def step_4_analysis(file_names, gff_file_name, genome_file):
    with open(gff_file_name) as f_in:
        with open(f"{gff_file_name[0:-3]}.rename.fa", "w") as f_out:
            for line in f_in:
                line = line.replace(">", ">Chr")
                f_out.writelines(line)
    cmd = f"gffread {gff_file_name[0:-3]}.rename.fa -g {genome_file} -x cds.fa"
    os.system(cmd)
    # -y protein.fa

    # gff_file_name = "TAIR10_GFF3_genes.gff"
    # extract_cds(genome_file, gff_file_name)
    cds_dict, chr_dict, strand_dict = get_cds_dictionary(gff_file_name)
    cds_sequence_dict = get_cds_sequence("cds.fa")

    codon_reads_dict = {}
    codon_counts_dict = {}
    with open(file_names, 'r') as f:
        for file_name in f:
            file_name = file_name.strip()
            bed_dic = get_bed_dict(file_name)

            for t_id in cds_sequence_dict:

                if not t_id.strip().endswith('.1'):
                    continue

                cds_sequence = cds_sequence_dict[t_id]
                if len(cds_sequence) % 3 != 0:
                    continue

                for i in range(0, len(cds_sequence) // 3):
                    codon = cds_sequence[int(i * 3):int(i * 3 + 3)]
                    transcript_location = i * 3

                    if codon not in codon_counts_dict:
                        codon_counts_dict[codon] = 1
                    else:
                        codon_counts_dict[codon] += 1

                    cds_list = cds_dict[t_id]
                    strand = strand_dict[t_id]
                    chr = chr_dict[t_id].strip("Chr")

                    transcript_location = transcript_location + 1
                    transcript_location = transcript_location - 16

                    chr_location = transcript_location_to_chr_location(transcript_location, strand, cds_list)

                    if strand == "+":
                        bed_key_1 = "%s_%s_%s" % (chr, int(chr_location + 0), strand)
                        bed_key_2 = "%s_%s_%s" % (chr, int(chr_location + 1), strand)
                        bed_key_3 = "%s_%s_%s" % (chr, int(chr_location + 2), strand)
                    else:
                        bed_key_1 = "%s_%s_%s" % (chr, int(chr_location - 0), strand)
                        bed_key_2 = "%s_%s_%s" % (chr, int(chr_location - 1), strand)
                        bed_key_3 = "%s_%s_%s" % (chr, int(chr_location - 2), strand)

                    reads = bed_dic.get(bed_key_1, 0) + bed_dic.get(bed_key_2, 0) + bed_dic.get(bed_key_3, 0)
                    reads = min(50, reads)

                    if codon not in codon_reads_dict:
                        codon_reads_dict[codon] = reads
                    else:
                        codon_reads_dict[codon] += reads

            # plot_counts_count(codon_counts_dict, file_name)
            # plot_reads_count(codon_reads_dict, file_name)
            plot_norm_count(codon_reads_dict, codon_counts_dict, file_name)
