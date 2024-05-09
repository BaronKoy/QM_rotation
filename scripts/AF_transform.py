complete_file = open('genstart.txt', 'w')
print('Position,AF', file = complete_file)
with open('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/SLiM_AFs/genstart.txt', 'r') as firstfile:
    for line in firstfile.readlines():
        if line.startswith(tuple('0123456789')):
            line = line.strip()
            line = line.split(sep=None, maxsplit=-1)
            locat = line[3]
            af1 = float(line[8])
            af2 = (af1 / 100)
            print(f'{locat},{af2}', file = complete_file)
firstfile.close()