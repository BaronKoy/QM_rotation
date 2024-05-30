# Convert Variant Call Format file to file containing location & minor allele frequency for plotting
import json

final = open('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/vcfs/cage_1/afs/gen2.csv', 'w')
print(f'Location,Minor_allele_frequency', file = final)
with open('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/vcfs/cage_1/gen2.vcf') as r:
    for line in r:
        if not line.startswith('#'):
            parts = line.split()
            info = parts[7]; info = info.split(';')
            ac = info[6]; ac = ac.split('='); ac_num = int(ac[1])
            freq = (ac_num/96)
            print(parts[1], freq, sep = ',', file = final)
r.close()