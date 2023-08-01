# launch tb-profiler env
import argparse
import os
import sys
import subprocess
import string
import random
import csv
import pandas as pd
import json
from Bio import SeqIO
import time


def file_choices(choices, fname):
    ext = os.path.splitext(fname)[1][1:]
    if ext not in choices:
        parser.error("file doesn't end with one of {}".format(choices))
    return fname


pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser(description="Prediction from csv file.")
parser.add_argument("-i", "--input", metavar="fastainput", help="Input fasta file(s)",
                    type=lambda s: file_choices(("fasta", "fna"), s), nargs='+')
parser.add_argument("-o", "--output", help="Directs the output to a name of your choice", type=str)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.1')


args = parser.parse_args()
in_file = args.input
output_file = args.output

if not os.path.exists("sitvit_geno/"):
    os.makedirs("sitvit_geno/")
else:
    for filename in os.listdir("sitvit_geno/"):
        os.remove("sitvit_geno/" + filename)

if os.path.exists("./sitvit_geno_final.tsv"):
    os.remove("./sitvit_geno_final.tsv")
    g = open("./sitvit_geno_final.tsv", 'w')
    g.write('')
    g.write('FilePath\tFlags\tRunName\tSpoligoType(MiruHero)\tMiruType\tLineage(MiruHero)\n')
    g.close()
else:
    g = open("./sitvit_geno_final.tsv", 'w')
    g.write('')
    g.write('FilePath\tFlags\tRunName\tSpoligoType(MiruHero)\tMiruType\tLineage(MiruHero)\n')
    g.close()

bashfile = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
bashfile = '/tmp/' + bashfile + '.sh'

miru_Order = ['miru02', 'miru04', 'miru10', 'miru16', 'miru20', 'miru23', 'miru24', 'miru26', 'miru27', 'miru31',
              'miru39', 'miru40', 'Mtub04', 'ETRC', 'Mtub21', 'Qub11b', 'ETRA', 'Mtub29', 'Mtub30', 'ETRB', 'Mtub34',
              'Mtub39', 'QUB26', 'QUB4156']

ancien = 'Mtub04', 'ETRC', 'Mtub21', 'Qub11b', 'ETRA', 'Mtub29', 'Mtub30', 'ETRB', 'Mtub34', 'Mtub39', 'QUB26', 'QUB4156'
new = 'ETR-A', 'ETRB', 'ETR-C', 'QUB-11b', 'QUB-26', 'QUB-4156', 'Mtub04', 'Mtub21', 'Mtub29', 'Mtub30', 'Mtub34', 'Mtub39'

# folder path
dir_path_sum_tsv = r'sitvit_geno/'

dir_path_fasta = r'./'

# list to store files
res_sum_tsv = []
res_fasta = []
lineages = []
resistance = []
spol_mirureader = []

debut = time.perf_counter()
# Iterate directory
for f in args.input:
    res_fasta.append(f)
print(res_fasta)

for i in res_fasta:
    subprocess.call(['MiruHero', '-m24', i, "-o", "./sitvit_geno"])

for file in os.listdir(dir_path_sum_tsv):
    # check only text files
    if file.endswith('.summary.tsv'):
        res_sum_tsv.append(file)
print(res_sum_tsv)

for i in range(len(res_sum_tsv)):
    with open("./sitvit_geno/" + res_sum_tsv[i] + "", "r+") as myfile:
        content = myfile.readlines()
        with open("./sitvit_geno_final.tsv", "a") as outfile:
            if len(content) == 2:
                outfile.write(content[1])
            else:
                continue
    outfile.close()

pd.set_option('display.max_columns', None)

df = pd.read_table('sitvit_geno_final.tsv', index_col=False, dtype={'SpoligoType(MiruHero)': 'string'},
                   usecols=['FilePath', 'RunName', 'SpoligoType(MiruHero)', 'MiruType', 'Lineage(MiruHero)'])

df_spollineage = df[['RunName', 'SpoligoType(MiruHero)', 'MiruType']].copy()

print('df_spol', df_spollineage)
df_spollineage.to_csv("df_spollineage_entry.csv", sep=';', index=False, header= False)

os.system('java -jar SpolLineages/spollineages.jar -i df_spollineage_entry.csv -o out-df_spol.csv')

if os.path.exists("spo_gca.out"):
    os.remove("spo_gca.out")

g = open("./spo_gca.out", 'w')
g.write('')
g.write('FilePath\tSpoligoType\tSpoligoType(Spotyping)\n')
g.close()

for i in range(len(df)):
    print(res_fasta[i])
    # check only text files
    os.system('python3 MIRUReader/MIRUReader.py -p mirus -r ' + df.loc[i, 'FilePath'] + '> miru.txt')
    os.system('tb-profiler profile -f' + df.loc[i, 'FilePath'] + '')
    os.system('python3 SpoTyping/SpoTyping-v3.0-commandLine/SpoTyping.py --noQuery --seq ' + df.loc[
        i, 'FilePath'] + ' -o spo_gca.out')

    df3 = pd.read_table('miru.txt', index_col=False)

    spol_mirureader.append(df3.iloc[0][1:].to_numpy())

    with open("./results/tbprofiler.results.json", "r") as f:
        data = json.load(f)

        if len(data['sublin']) == 0:
            lineages.append('ND')
        else:
            lineages.append(data['sublin'])

        if len(data['dr_variants']) == 0:
            resistance.append('ND')

        elif len(data['dr_variants']) > 0 and data['dr_variants'][0]['drugs'][0]['confers'] == 'sensitive':
            resistance.append("" + data['dr_variants'][0]['drugs'][0]['drug'] + "(S)")

        elif len(data['dr_variants']) > 0 and data['dr_variants'][0]['drugs'][0]['confers'] == 'resistance':
            resistance.append("" + data['dr_variants'][0]['drugs'][0]['drug'] + "(R)")

        elif len(data['dr_variants']) > 0 and data['dr_variants'][0]['drugs'][0]['confers'] == 'other':
            resistance.append("" + data['dr_variants'][0]['drugs'][0]['drug'] + "(O)")

        elif len(data['dr_variants']) > 0 and len(data['drtype']) == 0:
            resistance.append(data['dr_variants'][0]['drugs'][0]['drug'])

        f.close()
    os.remove('./results/tbprofiler.results.json')

df_csv_spol = pd.read_table('out-df_spol.csv', index_col=False, sep=';')
print('df csv spol', df_csv_spol["SubLineage (binary rules)"])

df2 = pd.read_table('spo_gca.out', index_col=False, dtype={'SpoligoType(Spotyping)': 'string'}, )
df.insert(loc=3, column='SpoligoType(Spotyping)', value=df2["SpoligoType(Spotyping)"])
df = df.join(df_csv_spol["SubLineage (binary rules)"])
# df.insert(loc=6, column='Family(SpolLineages)', value=df_csv_spol["SubLineage (binary rules)"])
df['Lineages'] = lineages
df['Resistance'] = resistance

for i in range(len(df)):
    mylist = df['MiruType'][i]
    myorder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 19, 13, 15, 22, 23, 12, 14, 17, 18, 20, 21]
    orderMIRUREader = [1, 4, 6, 7, 9, 15, 16, 17, 18, 20, 24, 5, 11, 14, 3, 10, 22, 23, 2, 8, 12, 13, 19, 21]

    mylist = [mylist[i] for i in myorder]
    s = ''.join(mylist)
    df['MiruType'][i] = s

for x in range(len(spol_mirureader)):
    for y in range(len(spol_mirureader[x])):
        if spol_mirureader[x][y] == 'ND':
            spol_mirureader[x][y] = '-'
        elif spol_mirureader[x][y] == 10:
            spol_mirureader[x][y] = 'A'
        elif spol_mirureader[x][y] == 11:
            spol_mirureader[x][y] = 'B'
        elif spol_mirureader[x][y] == 12:
            spol_mirureader[x][y] = 'C'
        elif spol_mirureader[x][y] == 13:
            spol_mirureader[x][y] = 'D'
        elif spol_mirureader[x][y] == 14:
            spol_mirureader[x][y] = 'E'
        elif spol_mirureader[x][y] == 15:
            spol_mirureader[x][y] = 'F'
        elif str(spol_mirureader[x][y]).count('s') > 0:
            index = str(spol_mirureader[x][y]).find('s')
            spol_mirureader[x][y] = spol_mirureader[x][y][:index]

new_miru = []
for i in range(len(spol_mirureader)):
    new_miru.append(''.join([str(v) for v in spol_mirureader[i]]))

values = ''.join([str(v) for v in spol_mirureader])
values2 = values.replace("'", '')
values2 = values2.replace(" ", "")

df["MiruType(mirureader)"] = pd.Series(new_miru)
df.rename(columns={'MiruType': 'MiruType(MiruHero)', 'Resistance': 'Resistance(tb-profiler)',
                   'Lineages': 'Lineages(tb-profiler)', "SubLineage (binary rules)":'Family(SpolLineages)'}, inplace=True)
df = df[["FilePath", "RunName", "SpoligoType(MiruHero)", "SpoligoType(Spotyping)", "MiruType(MiruHero)",
         "MiruType(mirureader)", "Lineage(MiruHero)", "Lineages(tb-profiler)", "Resistance(tb-profiler)", "Family(SpolLineages)"
         ]]

print(df)
df.to_csv("" + output_file + ".csv", index=False)
fin = time.perf_counter()
print(f" The runtime is {fin - debut:0.4f} seconds")
subprocess.call(['rm', '-rf', "*.fasta.*"])
