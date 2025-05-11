# Autor Maciej Szyszko
# Skrypt do pre-processingu - przerabianie pliku .struct w odpowiedni format do obliczeń struktur mieszanych (MoWS2) w Wien2K

import re 

ISPLIT = 8
# input_file -> nazwa pliku wejściowego
input_file = 'MoWS2.struct'
output_file = input_file + "_modified"

W_INDEX = 1

def edit_part(output, lines):
    global W_INDEX
    if lines[-4][0] == 'W':
        if len(lines) == 6:
            lines[2] = f'W{W_INDEX}'.ljust(10) + lines[2][10:]
            W_INDEX += 1
            output.writelines(lines)               
        else:
            ll = lines[:2] + lines[-4:]
            ll[1] = re.sub(r'MULT=\s*\d+', lambda m: 'MULT= 1', ll[1])
            ll[2] = f'W{W_INDEX}'.ljust(10) + ll[2][10:]
            W_INDEX += 1
            output.writelines(ll)
            for i in range(len(lines)-6):
                ll[0] = 'ATOM  ' + lines[2+i].lstrip()
                ll[2] = f'W{W_INDEX}'.ljust(10) + ll[2][10:]
                W_INDEX += 1
                output.writelines(ll)
    else:
        output.writelines(lines)

with open(input_file, 'r') as input:
    lines = input.readlines()
    with open(output_file, 'w') as output:
        for i in range(4):
            output.write(lines[i])

        start_i = 0
        end_i = 0
        for i in range(4, len(lines)):            
            if lines[i][:4] == 'ATOM':
                start_i = i
            if start_i != 0 and lines[i][:5] == 'LOCAL':
                end_i = i+2
            if start_i != 0 and end_i != 0:
                edit_part(output, lines[start_i:end_i+1])
                start_i = 0
                end_i = 0
            
            if 'NUMBER OF SYMMETRY OPERATIONS' in lines[i]:
                #jeśli nie chcesz kopiować operacji symetrii to zakomentuj linię poniżej (51)
                edit_part(output, lines[i:len(lines)])
                break
            

with open(output_file, 'r+') as output:
    edit_output = output.readlines()
    atom_n = 0
    for i in range(len(edit_output)):
        if edit_output[i][:4] == 'ATOM':
            atom_n += 1
            edit_output[i] = f'ATOM' + str(-atom_n).rjust(4) + edit_output[i][8:]
        elif edit_output[i][8] == ':':
            edit_output[i] = str(-atom_n).rjust(8) + edit_output[i][8:]
        elif 'ISPLIT' in edit_output[i]:
            edit_output[i] = re.sub(r'ISPLIT=\s*\d+', lambda m: 'ISPLIT= ' + str(ISPLIT), edit_output[i])
    edit_output[1] = re.sub(r'NONEQUIV.ATOMS\s*\d+', lambda m: 'NONEQUIV.ATOMS ' + str(atom_n).rjust(3), edit_output[1])
    output.seek(0)
    output.writelines(edit_output)
