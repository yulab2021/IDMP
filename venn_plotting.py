import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import csv

s1 = set()
s2 = set()
s3 = set()

with open("input.list", 'r') as f:
    for line in f:
        file1 = line.strip().strip("Aligned.sortedByCoord.out.bed") + "_co_translation.csv"
        file2 = line.strip().strip("Aligned.sortedByCoord.out.bed") + "_EJC_upstream_50.csv"
        file3 = line.strip().strip("Aligned.sortedByCoord.out.bed") + "uorf.csv"
        with open(file1, 'r') as csvfile1:
            reader1 = csv.reader(csvfile1)
            column = [row[0].strip().split('.')[0]+".1" for row in reader1]
            s1 = s1.union(set(column))
        with open(file2, 'r') as csvfile2:
            reader2 = csv.reader(csvfile2)
            column = [row[0].strip().split('.')[0]+".1" for row in reader2]
            s2 = s2.union(set(column))
        with open(file3, 'r') as csvfile3:
            reader3 = csv.reader(csvfile3)
            column = [row[0].strip().split('.')[0]+".1" for row in reader3]
            s3 = s3.union(set(column))
data = [s1, s2, s3]

g = venn3(data, set_labels=('CTRD', 'EJC', "uORF"),  # 设置组名
          set_colors=("#098154", "#c72e29", '#563F77'),  # 设置圈的颜色，中间颜色不能修改
          alpha=0.6,  # 透明度
          normalize_to=1.0)
plt.savefig("venn.pdf")
