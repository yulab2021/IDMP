import argparse
from mode1 import *
from mode2 import *
from mode3 import *
from mode4 import *
from mode5 import *

if __name__ == "__main__":
    # Adding Command Line Parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="List file includes the names of bed file")
    parser.add_argument('--mirna_name', help='The name of the miRNA that you want to analysis')
    parser.add_argument('--gff3_file', help='Pre-miRNA and mature miRNA annotation file')
    parser.add_argument('--gff', help='The general feature format file')
    parser.add_argument('--genome', help='The genome .fa file')
    parser.add_argument('--mirna_fa', help='MiRNA sequence files')
    parser.add_argument('--ps_file', help='Results of miRNA sequence alignment')
    parser.add_argument('--gene_name', help='Genename annotation file')
    parser.add_argument('--sample_info', help='Sample information annotation files')
    parser.add_argument('--rna_id', help='Three-letter abbreviation of the species to be studied')
    parser.add_argument('--mode', help='Input the mode that you want to analysis')

    args = parser.parse_args()

    if args.mode == '1':
        print("mode1:co-translational RNA decay analysis")
        mode_1_step_1(args.input, args.gff)
        get_csv_file(args.input, args.gff, args.gene_name)
        f_list = pattern_3(args.input, args.gff)
        pattern_3_plotting(f_list)
        step_4_analysis(args.input, args.gff, args.genome)
        print("The analysis was completed and the results have been deposited in the current working directory!")
        exit()
    elif args.mode == '3':
        print('mode3:Exon junction complex analysis')
        with open(args.input, 'r') as f:
            for file in f:
                file = file.strip()
                bed_dic = get_bed_dict(file)
                get_ejc_csv_file(file, bed_dic, args.gff, args.gene_name)
                get_ejc_heatmap(file)
        get_file_for_plotting(args.input, args.gff, args.sample_info)
        ejc_plotting("new_EJC_upstream_50.csv")
        print("The analysis was completed and the results have been deposited in the current working directory!")
        exit()
    elif args.mode == '4':
        print('mode4:5\'P reads distribution along pre-miRNA analysis')
        mode4_analysis(args.gff3_file, args.mirna_name, args.input)
        print("The analysis was completed and the results have been deposited in the current working directory!")
        exit()
    elif args.mode == '5-1':
        print('mode5-1:Identification of miRNA cleavage sites')
        rna_id = args.rna_id.lower()
        extract_cds(args.genome, args.gff)
        extract_mirna(args.mirna_fa, args.rna_id)
        print('Please download cds.fa and mirna.fa and submit to the provided website,and get the analysis result!')
        exit()
    elif args.mode == '5-2':
        print('mode5-2:Identification of miRNA cleavage sites')
        mode5_analysis(args.input, args.gff, args.ps_file)
        reads_plotting(args.mirna_name, args.input, args.gff)
        print("The analysis was completed and the results have been deposited in the current working directory!")
        exit()
    elif args.mode == '2':
        print('mode2:5\'UTR uORF identification')
        gene_name_dic = get_gene_name_dict(args.gene_name)
        uorf_analysis(args.input, args.gff, gene_name_dic, args.genome)
        with open(args.input, 'r') as f:
            for file in f:
                get_uorf_heatmap(file)
        print("The analysis was completed and the results have been deposited in the current working directory!")
        exit()
