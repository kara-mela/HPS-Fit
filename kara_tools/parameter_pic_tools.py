"""
- modules for creating diagrames: x, y chemical elements, color plot for lattice parameters
"""

import numpy as np
import itertools as it
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcls
from matplotlib.ticker import NullLocator, IndexLocator, LinearLocator, FixedLocator, MultipleLocator
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from kara_tools import chem_elem, TUBAF

def standardize(actual, maxi, mini=0.):
    """
    [0,1]-Normierung eines Datensatzes
    """
    return (actual-mini)/(maxi-mini)


def get_gtf(r_A , r_B, r_O = 132.):
    """
    gtf = Goldschmidt tolerance factor
    """
    return (r_A + r_O)/(np.sqrt(2.)*(r_B + r_O))

    
def similar_filter(values):
    """
    - in: values, out: values2
    - comparing lines of values with same chemical composition
    - reducing lines to max. 3 (due to limited possibilities of imaging)
    - choosing line with lattice parameters: max, min and one arbitrary intermediate
    """
    print "Filtering similar entries..."
    # "delete" similar entries
    mycount = 0     # number of equal compounds
    values2 = []
    for n in range(len(values)-1):
        if values[n][0:2] == values[n+1][0:2]: # if same composition in this and next line
            mycount += 1 # increment
        else: # else: take last block of same composition and manipulate
            dummy = values[n-mycount:n+1]
            
            pop_list = []       # line index of line that is similar to another one (can be popped)
            
            # choose entry-index k in dummy
            mylist = range(mycount+1)
            mylist2 = range(mycount+1)
            
            if len(mylist2) > 3: # more than three entries
                dummy_a = []
                dummy_super = []
                for entry in mylist2:
                    dummy_a.append(dummy[entry][3])     # list of a-parameters of independet entries
                    dummy_super.append(dummy[entry][6]) # list of superlattice values
                
                # determine max/min (and append them) and argmax/argmin
                a_max = max(dummy_a)
                a_argmax = dummy_a.index(a_max)
                values2.append(dummy[a_argmax])
                a_min = min(dummy_a)
                a_argmin = dummy_a.index(a_min)
                values2.append(dummy[a_argmin])
                
                dummy_index = []
                dummy_range = range(len(dummy_super))
                for n in range(len(dummy_super)):
                    if (dummy_super[n] != dummy_super[a_argmin]) and (dummy_super[n] != dummy_super[a_argmax]):
                        dummy_index.append(n)
                
                if len(dummy_index) > 0:
                    values2.append(dummy[dummy_index[0]])
                else:
                    for n in range(len(dummy_super)):
                        if (n != a_argmin) and (n != a_argmax):
                            dummy_index.append(n)
                    values2.append(dummy[dummy_index[0]])
            else:
                for entry in mylist2:
                    values2.append(dummy[entry])
            len(mylist2)
            
            mycount = 0
    
    return values2

    
def LaTeX_output(full, filtered, header, name):
    """
        - in: full data, filtered data, out: LaTeX-like document
    """
    print "Creating LaTeX-output..."
    # Ausgabe als LaTeX-freundliches *.txt
    mymarkers = []
    f = open(name, 'w')
    
    full = full.tolist()
    filtered = filtered.tolist()
    
    diff = lambda l1, l2: [x for x in l1 if x not in l2]
    marked = diff(full, filtered)
    
    k = 0
    for i in range(len(full)):
        for j in range(len(filtered)):
            if set(full[i]).issubset(set(filtered[j])):
            #if set(values[nr]).issubset(set(values2[nr2])):
            # if full[i] in filtered[j]:
                if 'unique' not in full[i][5]:
                    full[i][5] = 'unique' + full[i][5]
                    k += 1
        values_line = full[i]
        
        if 'unique' not in values_line[5]:
            for k in range(len(values_line)):
                values_line[k] = '\comment{' + str(values_line[k]) + '}'
        else:
            values_line[5] = values_line[5].replace('unique', '')
        
        for nr4 in range(len(values_line)):
            if header[nr4] == 'Autor':
                if '\comment' in values_line[nr4]:
                    values_line[nr4] = values_line[nr4].replace('\comment{', '\comment{\cite{')
                    values_line[nr4] += '}'
                else:
                    values_line[nr4] = '\cite{' + values_line[nr4] + '}'
                if '.' in values_line[nr4]:
                    values_line[nr4] = values_line[nr4].replace('.', ',')
                values_line[nr4] += ' & \t'
                
            elif header[nr4] == 'ICSD':
                values_line[nr4] += ' \\\\'
            
            elif header[nr4] == 'Kommentar':
                values_line[nr4] = ''
            
            elif header[nr4] == 'Ausdehnung':
                if len(values_line[nr4]) > 0:
                    if '1x1x1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('1x1x1', '$1 \\times 1 \\times 1$')
                    elif '2x2x2' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('2x2x2', '$2 \\times 2 \\times 2$')
                    elif '2x2x1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('2x2x1', '$2 \\times 2 \\times 1$')
                    elif '2x2x4' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('2x2x4', '$2 \\times 2 \\times 4$')
                    elif '2x2x?' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('2x2x?', '$2 \\times 2 \\times ?$')
                    elif '2x2x8' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('2x2x8', '$2 \\times 2 \\times 8$')
                    elif '?x?x1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('?x?x1', '$? \ \\times \ ? \ \\times 1$')
                    elif '???' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('???', 'N/A')
                    elif '6x2*sqrt(3)x4' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('6x2*sqrt(3)x4', '$6 \\times 2\\cdot \\sqrt{3} \\times 4$')
                    elif '6x2*sqrt(3)x2' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('6x2*sqrt(3)x2', '$6 \\times 2\\cdot\\sqrt{3} \\times 2$')
                    elif '1xsqrt(3)x1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('1xsqrt(3)x1', '$1 \\times \\sqrt{3} \\times 1$')
                values_line[nr4] += ' & \t'
            
            elif header[nr4] == 'Struktur':
                if len(values_line[nr4]) > 0:
                    if ',' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace(',', '')
                    #############
                    if values_line[nr4][0] == ' ':
                        values_line[nr4] = values_line[nr4][1:]
                    #############
                    if 'b = ' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('b = ', '$b = \ $')
                    #############
                    if 'P -1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P -1', '$P \ \overline{1} \ (2)$')
                        mymarkers.append('a')
                    elif 'P 1 2/m 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 1 2/m 1', '$P \ 1 \ 2/m \ 1 \ (10)$')
                        mymarkers.append('m')
                    elif 'P 1 21/m 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 1 21/m 1', '$P \ 1 \ 2_1/m \ 1 \ (11)$')
                        mymarkers.append('m')
                    elif 'C 1 2/m 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('C 1 2/m 1', '$C \ 1 \ 2/m \ 1 \ (12)$')
                        mymarkers.append('m')
                    elif 'P 1 21/c 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 1 21/c 1', '$P \ 1 \ 2_1/c \ 1 \ (14)$')
                        mymarkers.append('m')
                    elif 'C 1 2/c 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('C 1 2/c 1', '$C \ 1 \ 2/c \ 1 \ (15)$')
                        mymarkers.append('m')
                    elif 'orthorhombic' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('orthorhombic', '$o$')
                        mymarkers.append('o')
                    elif 'P 2 2 21' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 2 2 21', '$P \ 2 \ 2 \ 2_1 \ (17)$')
                        mymarkers.append('o')
                    elif 'P 21 21 21' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 21 21 21', '$P \ 2_1 \ 2_1 \ 2_1 \ (19)$')
                        mymarkers.append('o')
                    elif 'P m c 21' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P m c 21', '$P \ m \ c \ 2_1 \ (26)$')
                        mymarkers.append('o')
                    elif 'P c a 21' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P c a 21', '$P \ c \ a \ 2_1 \ (29)$')
                        mymarkers.append('o')
                    elif 'P b a 2' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P b a 2', '$P \ b \ a \ 2 \ (32)$')
                        mymarkers.append('o')
                    elif 'P n a 21' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P n a 21', '$P \ n \ a \ 2_1 \ (33)$')
                        mymarkers.append('o')
                    elif 'A m m 2' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('A m m 2', '$A \ m \ m \ 2 \ (38)$')
                        mymarkers.append('o')
                    elif 'P m m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P m m m', '$P \ m \ m \ m \ (47)$')
                        mymarkers.append('o')
                    elif 'P n n a' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P n n a', '$P \ n \ n \ a \ (52)$')
                        mymarkers.append('o')
                    elif 'P b a m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P b a m', '$P \ b \ a \ m \ (55)$')
                        mymarkers.append('o')
                    elif 'P b c m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P b c m', '$P \ b \ c \ m \ (57)$')
                        mymarkers.append('o')
                    elif 'P m m n Z' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P m m n Z', '$P \ m \ m \ n \ Z \ (59)$')
                        mymarkers.append('o')
                    elif 'P n m a' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P n m a', '$P \ n \ m \ a \ (62)$')
                        mymarkers.append('o')
                    elif 'C m c m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('C m c m', '$C \ m \ c \ m \ (63)$')
                        mymarkers.append('o')
                    elif 'C m m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('C m m m', '$C \ m \ m \ m \ (65)$')
                        mymarkers.append('o')
                    elif 'F m m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('F m m m', '$F \ m \ m \ m \ (69)$')
                        mymarkers.append('o')
                    elif 'F d d d' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('F d d d', '$F \ d \ d \ d \ (70)$')
                        mymarkers.append('o')
                    elif 'I m m a' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('I m m a', '$I \ m \ m \ a \ (74)$')
                        mymarkers.append('o')
                    elif 'tetragonal' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('tetragonal', '$t$')
                        mymarkers.append('t')
                    elif 'I 4/m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('I 4/m', '$I \ 4/m \ (87)$')
                        mymarkers.append('t')
                    elif 'P 4 m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 4 m m', '$P \ 4 \ m \ m \ (99)$')
                        mymarkers.append('t')
                    elif 'P 4/m b m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 4/m b m', '$P \ 4/m \ b \ m \ (127)$')
                        mymarkers.append('t')
                    elif 'P 4/m m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 4/m m m', '$P \ 4/m \ m \ m \ (139)$')
                        mymarkers.append('t')
                    elif 'I 4/m c m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('I 4/m c m', '$I \ 4/m \ c \ m \ (140)$')
                        mymarkers.append('t')
                    elif 'I 41/a m d' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('I 41/a m d', '$I \ 4_1/a \ m \ d \ (141)$')
                        mymarkers.append('t')
                    elif 'R -3 H' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('R -3 H', '$R \ \overline{3} \ H \ (148)$')
                        mymarkers.append('3')
                    elif 'R 3 m H' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('R 3 m H', '$R \ 3 \ m \ H \ (160)$')
                        mymarkers.append('3')
                    elif 'R 3 c H' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('R 3 c H', '$R \ 3 \ c \ H \ (161)$')
                        mymarkers.append('3')
                    elif 'P -3 m 1' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P -3 m 1', '$P \ \overline{3} \ m \ 1 \ (164)$')
                        mymarkers.append('3')
                    elif 'R -3 m H' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('R -3 m H', '$R \ \overline{3} \ m \ H \ (166)$')
                        mymarkers.append('3')
                    elif 'R -3 c H' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('R -3 c H', '$R \ \overline{3} \ c \ H \ (167)$')
                        mymarkers.append('3')
                    elif 'hexagonal' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('hexagonal', '$h$')
                        mymarkers.append('h')
                    elif 'Ce2CoSi3' in values_line[nr4]:
                        # values_line[nr4] = values_line[nr4].replace('Ce2CoSi3', '\ce{Ce2CoSi3}')
                        values_line[nr4] = values_line[nr4].replace('Ce2CoSi3', '$P \ 6/m \ m \ m \ (191)$')
                        mymarkers.append('h')
                    elif 'P 63 c m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 63 c m', '$P \ 6_3 \ c \ m \ (185)$')
                        mymarkers.append('h')
                    elif 'P -6 2 c' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P -6 2 c', '$P \ \overline{6} \ 2 \ c \ (190)$')
                        mymarkers.append('h')
                    elif 'P 6/m m m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 6/m m m', '$P \ 6/m \ m \ m \ (191)$')
                        mymarkers.append('h')
                    elif 'P 63/m m c' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P 63/m m c', '$P \ 6_3/m \ m \ c \ (194)$')
                        mymarkers.append('h')
                    elif 'P n -3 Z' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P n -3 Z', '$P \ n \ \overline{3} \ Z \ (201)$')
                        mymarkers.append('c')
                    elif 'I a -3' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('I a -3', '$I \ a \ \overline{3} \ (206)$')
                        mymarkers.append('c')
                    elif 'P m -3 m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('P m -3 m', '$P \ m \ \overline{3} \ m \ (221)$')
                        mymarkers.append('c')
                    elif 'F m -3 m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('F m -3 m', '$F \ m \ \overline{3} \ m \ (225)$')
                        mymarkers.append('c')
                    elif 'F d -3 m Z' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('F d -3 m Z', '$F \ d \ \overline{3} \ m \ Z \ (227)$')
                        mymarkers.append('c')
                    elif 'F d -3 m' in values_line[nr4]:
                        values_line[nr4] = values_line[nr4].replace('F d -3 m', '$F \ d \ \overline{3} \ m \ Z \ (227)$')
                        mymarkers.append('c')
                    else:
                        mymarkers.append('')
                values_line[nr4] += ' & \t'
            
            
            else:
                values_line[nr4] += ' & \t'
            
            f.write(values_line[nr4])
        f.write('\n')
    
    f.close()
    
    mymarkers = np.array(mymarkers)
    
    return mymarkers

    
def decrypt_elements(values2, substance):
    """
    - in: filtered indata, atomic number list, out: array of lattice parameters, atomic numbers
    - replacing element symbol by atomic number (iterable) for A -> xs, B -> ys
    - creating lists for lattice parameters
    """
    a, b, c, xs, ys, gtf, ca, super  = [], [], [], [], [], [], [], []
    dummy_hexa, dummy_tetra, dummy_hexa_ca, dummy_tetra_ca = [], [], [], []
    
    elements = chem_elem.elem_list()
    
    # atomic numbers of I: (A2+)(B4+)O3, II: (A3+)(B3+)O3, III: (A1+)(B5+)O3
    A_I   = (12, 20, 38, 56, 82)
    A_II  = range(57, 72)
    A_III = (3, 11, 19, 37)
    B_I   = (22, 40, 72, 50, 58)
    B_II  = (13, 21, 23, 24, 25, 26, 27)
    B_III = (41, 73, 51)
    
    for row in values2:
        # attaching the chemical symbol to atomic number for A and B/ R and T...
        for Z in range(len(elements)):
            # element A resp. R
            if row[0] == elements[Z][0]:
                xs.append(Z)
            # element B resp. T
            if row[1] == elements[Z][0]:
                ys.append(Z)
                
        if substance == 'STO':
            # finding oxidation states and ionic radii
            if xs[-1] in A_I:
                ox_A = '2'
            elif xs[-1] in A_II:
                ox_A = '3'
            elif xs[-1] in A_III:
                ox_A = '1'
                
            if ys[-1] in B_I:
                ox_B = '4'
            elif ys[-1] in B_II:
                ox_B = '3'
            elif ys[-1] in B_III:
                ox_B = '5'
            
            r_A = elements[xs[-1]][1][ox_A]
            r_B = elements[ys[-1]][1][ox_B]
            
            ratio = get_gtf(r_A, r_B)
            gtf.append(ratio)
        else:
            gtf = 0
            
        """
        relation to header!!!
        """
        a.append(float(row[3].replace('\comment{', '').replace('};', '')))
        b.append(float(row[4].replace('\comment{', '').replace('};', '')))
        c.append(float(row[5].replace('\comment{', '').replace('};', '')))
        ca.append(float(row[5].replace('\comment{', '').replace('};', '')))
        
        if '2x2x' in row[6]:
            row[6] = row[6].replace('2x2x', '')
            super.append(float(row[6]))
        elif '1x1x' in row[6]:
            super.append(0)
        else:
            super.append(-1)
        
        
        # differentiate between hexagonal and tetragonal structures
        if c[-1]/a[-1] > 2:
            dummy_tetra.append(c[-1])
            dummy_tetra_ca.append(ca[-1])
        else:
            dummy_hexa.append(c[-1])
            dummy_hexa_ca.append(ca[-1])
    
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    ca = np.array(ca)
    xs = np.array(xs)
    ys = np.array(ys)
    super = np.array(super)
    gtf = np.array(gtf)
    dummy_hexa = np.array(dummy_hexa)
    dummy_tetra = np.array(dummy_tetra)
    dummy_hexa = np.array(dummy_hexa_ca)
    dummy_tetra = np.array(dummy_tetra_ca)
    
    return a, b, c, xs, ys, gtf, ca, super, dummy_hexa, dummy_tetra, dummy_hexa_ca, dummy_tetra_ca

    
def arrange_data(path, criteria, header, substance):
    """
    - in: path to data, criteria[0] = necessary element (Si, O), 
      criteria[1] = lattice parameter restriction 
    - out: arrays for lattice parameters, atomic numbers, marker style
    - reading hole data set
    - removing empty header, empty lines, and incomplete data
    - sorting: alphabetical in A, alphabetical in B
    - conformizing space groups (equalizinig different settings)
    - filtering similar entries
    - producing LaTeX-output
    """
    print "Arranging data..."
    # read in
    f = open(path, 'r')
    data = f.read()
    f.close()

    # create 2d-array
    data2 = []
    data = data.rsplit('\n')            # line by line
    for line in data:
        line = line.replace(',', '.')   # numbers to english setting
        line2 = line.rsplit(';')        # column by column
        data2.append(line2)
    data2 = np.array(data2)
    
    # cutout unimportant lines (missing parameter, header, etc)
    values = []
    for line in data2[1:]:          # remove header
        if len(line) > 1:           # ignore empty lines
            if line[3] != '' and line[2] == criteria[0] and float(line[criteria[1][0]]) < criteria[1][1]:       # remove lines with non-oxides, missing parameters
                if line[1] == '':
                    line[1] = 'Si'
                values.append(line)

    # sort via elements (resultat: alphabetical in A, within this alphabetical in B
    values = sorted(values, key=lambda B: B[1])
    values = sorted(values, key=lambda A: A[0])
    
    values2 = similar_filter(values)
    
    values2 = np.array(values2)
    values = np.array(values)
    
    if substance == 'HPS':
        file_name = 'HPS-Parameter-LaTeX.tex'
    elif substance == 'STO':
        file_name == 'Perowskit-Parameter-LaTeX.tex'
    elif substance == 'YMFO':
        file_name == 'YMFO-Parameter-LaTeX.tex'
    
    mymarkers = LaTeX_output(values, values2, header, name=file_name)
    
    return values2, mymarkers

     
def set_symbols(param, data, xs, ys, mymarkers, ax, owncolor, tetra=0., hexa=0.):
    """
    - setting symbols for lattice parameters
    - different shape, depending on crystal class
    - different filling, depending on amount of different data
    """
    
    
    """
    combine with values and values2!!!
    """
    print "Setting symbols..."
    
    # filling in symbols
    val = []
    repeats = 0
    print len(data)
    for n in range(len(data)):
        if n%10 == 0: print n+1, "/", len(data)
        # setting fillstyle: only one dataset -> full, else iterating half-fill (left, right, bottom (top))
        if (n < len(data)-2):
            if (xs[n] == xs[n+1] and ys[n] == ys[n+1]):
                mystyle = mlines.Line2D.fillStyles[repeats+1]
                repeats += 1
            else:
                if repeats != 0:
                    mystyle = mlines.Line2D.fillStyles[repeats+1]
                else:
                    mystyle = mlines.Line2D.fillStyles[repeats]
                repeats = 0 
        else:
            if (xs[n] == xs[n-1] and ys[n] == ys[n-1]):
                mystyle = mlines.Line2D.fillStyles[repeats+1]
            else:
                mystyle = mlines.Line2D.fillStyles[0]
                
        for axis in ax.itervalues():
            # attaching markers
            if param == 'super':
                mymarker = 'o'
                if data[n] == 0:
                    mysize = 4
                elif data[n] == -1:
                    mymarker = 'x'
                    mysize = 8
                else:
                    mysize = 13
            else:
                if mymarkers[n] == 'o': # orthorhombic
                    mymarker = 'D'
                    mysize=7
                elif mymarkers[n] == 't': # tetragonal
                    mymarker = '|'
                    mysize=7
                elif mymarkers[n] == 'c': # cubic
                    mymarker = 's'
                    mysize=7
                elif mymarkers[n] == 'm' or mymarkers[n] == 'a': # monoclinic/triclinic
                    mymarker = 'o'
                    mysize=7
                elif mymarkers[n] == 'h': # hexagonal
                    mymarker = 'h'
                    mysize=10
                elif mymarkers[n] == '3': # trigonal
                    mymarker = '^'
                    mysize=10
                
            # attaching colors
            if param == 'c':
                if data[n] > 10:
                    val.append(int((data[n]-tetr.min())/(tetr.max()-tetra.min())*255.))
                else:
                    val.append(int((data[n]-hexa.min())/(hexa.max()-hexa.min())*255.))
                mycolor = owncolor(val[-1])
            elif param == 'c/a':
                if data[n] > 2:
                    val.append(int((data[n]-tetra.min())/(tetra.max()-tetra.min())*255.))
                else:
                    val.append(int((data[n]-hexa.min())/(hexa.max()-hexa.min())*255.))
                mycolor = owncolor(val[-1])
            elif param == 'super':
                if data[n] < 1:
                    mycolor = 'black'
                else:
                    val.append(int((data[n]-1)/(8-1)*255))
                    mycolor = owncolor(val[-1])
            else:
                val.append(int((data[n]-data.min())/(data.max()-data.min())*255.))
                mycolor = owncolor(val[-1])
                
            
            axis.plot(xs[n], ys[n], c=mycolor, markersize=mysize, fillstyle=mystyle, marker=mymarker, alpha=1., linewidth=2)

    
def set_diagram_size(x_range, y_range):
    breite, hohe = [], []
    
    for x_entry in x_range:
        breite.append(x_entry[1]+1 - x_entry[0] + 1)
    breite.append(0.6)  # for colorbar
    
    for y_entry in y_range:
        hohe.append(y_entry[1] - y_entry[0] + 1)
    
    return breite, hohe

    
def set_axes(gs, hohe, breite, x_range, y_range):
    print "Setting axes..."
    # setting axes
    ax = {}
    k = 0
    for i,j in it.product(range(1,len(hohe)+1), range(1, len(breite)+1)):
        if i == 1:
            if j == 1:
                ax[str(i) + str(j)] = plt.subplot(gs[k])
                k += 1
            elif j == len(breite):
                k+=1
            else:
                ax[str(i) + str(j)] = plt.subplot(gs[k], sharey=ax[str(i) + "1"])
                k+=1
        else:
            if j == 1:
                ax[str(i) + str(j)] = plt.subplot(gs[k], sharex=ax["1" + str(j)])
                k+=1
            elif j == len(breite):
                k+=1
            else:
                ax[str(i) + str(j)] = plt.subplot(gs[k], sharex=ax["1" + str(j)], sharey=ax[str(i) + "1"])
                k+=1
    
    for axis in ax.values():
        axis.grid(False)
    
    # formatting outter lines
    for pos in ax.iterkeys():
        for spine in ax[pos].spines.keys():
            if pos[1] == "1" and spine == "left":
                pass
            elif pos[1] == str(len(breite)-1) and spine == "right":
                pass
            elif pos[0] == "1" and spine == "top":
                pass
            elif pos[0] == str(len(hohe)) and spine == "bottom":
                pass
            else:
                ax[pos].spines[spine].set_visible(False)
    
    # formatting tick positions, relative to subplots
    for pos in ax.iterkeys():
        if pos[1] == "1":
            ax[pos].yaxis.tick_left()
        elif pos[1] == str(len(breite)-1):
            ax[pos].yaxis.tick_right()
        else:
            ax[pos].yaxis.set_ticks_position('none')
        if pos[0] == "1":
            ax[pos].xaxis.tick_top()
        elif pos[0] == str(len(hohe)):
            ax[pos].xaxis.tick_bottom()
        else:
            ax[pos].xaxis.set_ticks_position('none')
        
        ax[pos].tick_params(labeltop='off')
        ax[pos].tick_params(labelright='off')
        if pos[1] != "1":
            ax[pos].tick_params(labelleft='off')
        if pos[0] != str(len(hohe)):
            ax[pos].tick_params(labelbottom='off')
    
    # formatting tick range
    for pos in ax.iterkeys(): 
        if pos[0] == "1" and pos[0] != len(breite)-1:
            ax[pos].xaxis.set_major_locator(FixedLocator(range(x_range[int(pos[1])-1,0]-1, 
              x_range[int(pos[1])-1,1]+1)))
        if pos[1] == "1":
            ax[pos].yaxis.set_major_locator(FixedLocator(range(y_range[int(pos[0])-1,0]-1, 
              y_range[int(pos[0])-1,1]+1)))
    
    
    elements = chem_elem.elem_list()
    
    # setting chemical symbol beside the atomic number
    for pos in ax.iterkeys():
        if pos[1] == "1":
            ax[pos].set_yticklabels(map(lambda y: "%s %i"%(elements[int(y)][0], int(y)), 
              ax[pos].get_yticks()))
        if pos[0] == "1" and pos[0] != len(breite)-1:
            ax[pos].set_xticklabels(map(lambda x: "%i \n %s"%(int(x), elements[int(x)][0]), 
              ax[pos].get_xticks()))
    
    # range of subplots
    for pos in ax.iterkeys():
        if pos[0] == "1":
            ax[pos].set_xlim(x_range[int(pos[1])-1,0] - 0.5, x_range[int(pos[1])-1,1] + 0.5)
        if pos[1] == "1" and pos[0] != len(breite)-1:
            ax[pos].set_ylim(y_range[int(pos[0])-1,0] - 0.5, y_range[int(pos[0])-1,1] + 0.5)
    
    # axis gap
    d = 0.17
    for pos in ax.iterkeys(): 
        if pos[0] == "1":
            kwargs = dict(transform=ax[pos].transAxes, color='k', clip_on=False)
            if pos[1] != str(len(breite)-1):
                ax[pos].plot((1-d/breite[int(pos[1])-1],1+d/breite[int(pos[1])-1]), 
                  (1-d/hohe[int(pos[0])-1],1+d/hohe[int(pos[0])-1]), **kwargs) # line at top right
            if pos[1] != "1":
                ax[pos].plot((-d/breite[int(pos[1])-1],d/breite[int(pos[1])-1]), 
                  (1-d/hohe[int(pos[0])-1],1+d/hohe[int(pos[0])-1]), **kwargs) # line at top left
        if pos[0] == str(len(hohe)):
            kwargs = dict(transform=ax[pos].transAxes, color='k', clip_on=False)
            if pos[1] != str(len(breite)-1):
                ax[pos].plot((1-d/breite[int(pos[1])-1],1+d/breite[int(pos[1])-1]), 
                  (-d/hohe[int(pos[0])-1],d/hohe[int(pos[0])-1]), **kwargs) # line at bottom right
            if pos[1] != "1":
                ax[pos].plot((-d/breite[int(pos[1])-1],d/breite[int(pos[1])-1]), 
                  (-d/hohe[int(pos[0])-1],d/hohe[int(pos[0])-1]), **kwargs) # line at bottom left
        
        if pos[1] == "1":
            kwargs = dict(transform=ax[pos].transAxes, color='k', clip_on=False)
            if pos[0] != str(len(hohe)):
                ax[pos].plot((-d/breite[int(pos[1])-1],d/breite[int(pos[1])-1]), 
                  (-d/hohe[int(pos[0])-1],d/hohe[int(pos[0])-1]), **kwargs) # line at bottom left
            if pos[0] != "1":
                ax[pos].plot((-d/breite[int(pos[1])-1],d/breite[int(pos[1])-1]), 
                  (1-d/hohe[int(pos[0])-1],1+d/hohe[int(pos[0])-1]), **kwargs) # line at top left

        if pos[1] == str(len(breite)-1):
            kwargs = dict(transform=ax[pos].transAxes, color='k', clip_on=False)
            if pos[0] != str(len(hohe)):
                ax[pos].plot((1-d/breite[int(pos[1])-1],1+d/breite[int(pos[1])-1]), 
                  (-d/hohe[int(pos[0])-1],d/hohe[int(pos[0])-1]), **kwargs) # line at bottom right
            if pos[0] != "1":
                ax[pos].plot((1-d/breite[int(pos[1])-1],1+d/breite[int(pos[1])-1]), 
                  (1-d/hohe[int(pos[0])-1],1+d/hohe[int(pos[0])-1]), **kwargs) # line at top right
    
    return ax

    
def plotting(param, data, xs, ys, mymarkers, my_style, substance, x_range, y_range, tetra=0., hexa=0.):
    print "Start plotting..."
    # color settings
    width = TUBAF.width(my_style)
    name = TUBAF.name(my_style)
    owncolor = mcls.LinearSegmentedColormap.from_list('own', 
       [TUBAF.rot(my_style), TUBAF.orange(my_style), TUBAF.gelb(my_style), 
       TUBAF.gruen(my_style), TUBAF.blau(my_style), TUBAF.cyan(my_style)])
    
    # size, gride ranges
    breite, hohe = set_diagram_size(x_range, y_range)
    
    f = plt.figure()
    gs = gridspec.GridSpec(len(hohe), len(breite), width_ratios=breite, height_ratios=hohe)
    
    ax = set_axes(gs, hohe, breite, x_range, y_range)
    
    set_symbols(param, data, xs, ys, mymarkers, ax, owncolor, tetra=0., hexa=0.)
    
    # # background color
    # ax["81"].set_axis_bgcolor("black") # A2+ B4+
    # # A3+ B3+ (lanthanides)
    # # A1+ B5+
    
    print "Setting legend..."
    # legend
    if param == 'a':
        mynorm = mcls.Normalize(vmin = min(a), vmax = max(a))
        mylabel = '$a$ [$\\AA$]'
    elif param == 'b':
        mynorm = mcls.Normalize(vmin = min(b), vmax = max(b))
        mylabel = '$b$ [$\\AA$]'
    elif param == 'c':
        mynorm = mcls.Normalize(vmin = min(c), vmax = max(c))
        mylabel = '$c$ [$\\AA$]'
    elif param == 'c/a':
        mynorm = mcls.Normalize(vmin = min(ca), vmax = max(ca))
        mylabel = 'ratio $c/a$'
    elif param == 'gtf':
        mynorm = mcls.Normalize(vmin = min(gtf), vmax = max(gtf))
        mylabel = 'Goldschmidt\'s tolerance factor $t$'
    elif param == 'super':
        mynorm = mcls.Normalize(vmin=1.01, vmax=8)
        mylabel = 'multiplicity in $c$'
    
    axim = plt.subplot(gs[:, -1])
    cb = matplotlib.colorbar.ColorbarBase(axim, cmap = owncolor, norm = mynorm)
    
    if param == 'super':
        cb.set_ticks(np.arange(0., 8., 1.))
        
    # distances between subplots
    print "Adjusting subplots.."
    plt.subplots_adjust(wspace=0.07, hspace=0.07, left=0.11, bottom=0.12, right=0.93, top=0.98)

    # labelling
    if substance == 'STO':
        f.text(0.5, 0.02, '$A$, chemical symbol', ha='center', va='center', fontsize=14)
        f.text(0.02, 0.5, '$B$, chemical symbol', ha='center', va='center', rotation='vertical', 
          fontsize=14)
    else:
        f.text(0.5, 0.02, '$R$, chemical symbol', ha='center', va='center', fontsize=14)
        f.text(0.02, 0.5, '$T$, chemical symbol', ha='center', va='center', rotation='vertical', 
          fontsize=14)
    
    f.text(0.98, 0.5, mylabel, ha='center', va='center', rotation='vertical', fontsize=14)

    plt.savefig(str(param).replace('/', '') + '-' + name + '.pdf', transparent=True)
    plt.show()

