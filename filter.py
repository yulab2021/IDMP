import subprocess
import os

# Extract all .fq file names and make a list

file_names = []
for _, _, files in os.walk('.', topdown=True):
    for name in files:
        if '.fq' in name and "trimed" not in name:
            file_names.append(name)

for file_name in file_names:
    i = 0
    output = file_name.strip('.fq') + '_trimed_filter.fq'
    with open(file_name, 'r') as fq:
        with open(output, 'w') as w:
            for line in fq:
                i = i+1
                if i % 4 == 1:
                    tmp1 = line
                if i % 4 == 2:
                    tmp2 = line
                    if len(line) < 18:
                        print_flag = 0
                    else:
                        print_flag = 1
                if i % 4 == 3:
                    tmp3 = line
                if i % 4 == 0:
                    tmp4 = line
                    if print_flag == 1:
                        w.write(tmp1)
                        w.write(tmp2)
                        w.write(tmp3)
                        w.write(tmp4)

    subprocess.run(["mv", f"{output}", "filtered_data"])
