import xml.dom.minidom
import re
import pandas as pd
import csv


def ConvertUnits(HFlist):
    for i in range(len(HFlist)):
        HFlist[i] = HFlist[i].upper()
        HFlist[i] = re.sub(',','',HFlist[i])
            
        if 'KCAL' in HFlist[i]:
            main_end = 0
            if '+' in HFlist[i]: # if uncertainty provided
                main_end = HFlist[i].find('+')
                un_begin = main_end + 3
                un_end = HFlist[i].find('K')
                uncertainty = HFlist[i][un_begin:un_end]
                uncertainty_convert = float(uncertainty) * 4.184
                convert_unc_part = '+/-' + str(uncertainty_convert)
            else: # if uncertatinty NaN
                main_end = HFlist[i].find('K')
                convert_unc_part = ''
            
            num_part = HFlist[i][:main_end]
            convert_num_part = float(num_part)*4.184
            
            HFlist[i] = str(round(convert_num_part,3)) + convert_unc_part + ' KJ'

        HFlist[i] = re.sub('NAN','NaN',HFlist[i])
        if ((HFlist[i] != 'NaN') and ('KJ' not in HFlist[i])):
            HFlist[i] = HFlist[i] + ' KJ'
    return HFlist


def ExtractMinMain(HFlist):
    mainvalues = []
    for i in range(len(HFlist)):
        unctt = HFlist[i].find('+')
        if unctt != -1 and HFlist[i] != 'NaN': # hf provided but no uncertainty  
            mainvalues.append(float(HFlist[i][:unctt]))
        elif unctt == -1 and HFlist[i] != 'NaN': # both hf and uncertainty provided
            findkj = HFlist[i].find('KJ')
            mainvalues.append(float(HFlist[i][:findkj]))
            
    minvalue = min(mainvalues)
    min_pos = mainvalues.index(minvalue)
    return [minvalue,min_pos]


def FindSubStr(substr, string, count): # return a list of all the positions
    index = 0
    result = []
    
    for j in range(count):
        next_substr_pos = string.find(substr,index)
        result.append(next_substr_pos)
        index = next_substr_pos+1
    return result


def Deconstruct(formula):
    if formula[-1] == ')': # XX(XX)
        formula = formula[:-1]
        last_b = formula.rfind('(')
        formula = formula[:last_b]+formula[last_b+1:]

    num_b = formula.count('(') 
    for bracket in range(num_b):
        b_pos1 = formula.find('(')
        b_pos2 = formula.find(')')
        
        front = formula[:b_pos1]
        inside = formula[b_pos1+1:b_pos2]
        
        if formula[b_pos2+1].isalpha(): # XX(XX)XX
            num = 1
            back = formula[b_pos2+1:]
        elif formula[b_pos2+1:].isdigit(): # XX(XX)n
            back = ''
            num = int(formula[b_pos2+1:])
        else: # XX(XX)nXX
            for i in range(1,10):
                if formula[b_pos2+i].isalpha():
                    back_start = b_pos2+i
                    break
            back = formula[b_pos2+i:]
            num = int(formula[b_pos2+1:b_pos2+i])
        
        formula = front
        for n in range(num):
            formula = formula+inside
        formula = formula+back
    return formula


def Addones(formula):
    l = list(formula)
    if l[-1] == '+' or l[-1] == '-': # If ion, the second last position must be digit. If not, add '1'
        if not l[-2].isdigit():
            l[-2] = l[-2] + '1'

    else:# If not ion, the last position must be digit. If not, add '1'
        if l[-1].isalpha():
            l[-1] = l[-1] + '1'

    for i in range(1,len(formula)): # add '1's after the elements which are not followed by digits
        if l[i].isupper() and not l[i-1].isdigit():
            l[i-1] = l[i-1] + '1'
            
    return ''.join(l)


def SumupElements(formula):
    elements = []
    final_formula = []
    formula_to_list = list(formula)
    for i in range(len(formula_to_list)):
        if formula_to_list[i].islower():
            formula_to_list[i-1] = formula_to_list[i-1] + formula_to_list[i]
            formula_to_list[i] = '' # move the seperate letters to form full elements
    for j in range(1,len(formula_to_list)):
        if formula_to_list[j].isdigit() and formula_to_list[j-1].isdigit():
            formula_to_list[j-1] = formula_to_list[j-1] + formula_to_list[j]
            formula_to_list[j] = ''    
            
    while '' in formula_to_list:
        formula_to_list.remove('')

    for e in formula_to_list:
        if (e.isalpha()) and (e not in elements): # e is a new element
            elements.append(e)
            elements.sort() # list of all the unique elements

    for element in elements:
        final_formula.append(element)
        total_num = 0
        for k in range(len(formula_to_list)):
            if formula_to_list[k] == element:
                total_num = total_num + int(formula_to_list[k+1])
        final_formula.append(str(total_num))
    if formula[-1] == '+' or formula[-1] == '-':
        final_formula.append(formula[-1])

    summed_formula = ''.join(final_formula)

    
    if summed_formula[-1] == '1' and not summed_formula[-2].isdigit():
        summed_formula = summed_formula[:-1]

    if len(summed_formula) > 2 and summed_formula[-2] == '1' and not summed_formula[-3].isdigit() and (summed_formula[-1] == '+' or summed_formula[-1] == '-'):
        summed_formula = summed_formula[:-2]+summed_formula[-1]
        
    formula_l = list(summed_formula)
    for k in range (1, len(formula_l)-1):
        if formula_l[k] == '1' and formula_l[k+1].isalpha() and formula_l[k-1].isalpha():
            formula_l[k]='#'
    while '#' in formula_l:
        formula_l.remove('#')
    result = ''.join(formula_l)
    return result




DOMTree = xml.dom.minidom.parse("BURCAT_THR.xml")
database = DOMTree.documentElement


species = database.getElementsByTagName("specie")

substance_l=[]
hflist=[]
highs=[]
lows=[]
ref_codes=[]
reflist=[]
reflist_adt=[]
hfvalues_a=[]

for specie in species:
    specie_high=[]
    specie_low=[]
    specie_ref_adt=[]
    
    if specie.getElementsByTagName('phase') != []:
        phases = specie.getElementsByTagName('phase')[0]
        status = phases.getElementsByTagName('phase')[0]
        
        if status.childNodes[0].data == 'G': # only need gas phase
            formula1 = specie.getElementsByTagName('formula_name_structure')
            formula2 = phases.getElementsByTagName('formula')[0]
            
            source = phases.getElementsByTagName('source')[0]

            hf298s = specie.getElementsByTagName('hf298')
            references = specie.getElementsByTagName('reference')
            additionals = specie.getElementsByTagName('additional_information')

            coeffs = phases.getElementsByTagName('coefficients')[0]

            high_coeffs = coeffs.getElementsByTagName('range_1000_to_Tmax')[0]
            high_coef = high_coeffs.getElementsByTagName('coef')

            low_coeffs = coeffs.getElementsByTagName('range_Tmin_to_1000')[0]
            low_coef = low_coeffs.getElementsByTagName('coef')
            
            
            if formula1.length != 0: # if <formula_name_structure> empty: take formula from <formula>
                for formul in formula1:
                    subs = formul.getElementsByTagName('formula_name_structure_1')[0]
                    sub_formula = subs.childNodes[0].data
                    substance_l.append(sub_formula)
            else:
                sub_formula = formula2.childNodes[0].data
                substance_l.append(sub_formula)

            ref_codes.append(source.childNodes[0].data)
            
            
            # references
            if references.length == 1: # if references provided
                for ref in references:
                    refs = []
                    ref1 = ref.getElementsByTagName('reference_1')[0]
                    refs.append(ref1.childNodes[0].data)
                    
                    if len(ref.getElementsByTagName('reference_2')) != 0:
                        ref2 = ref.getElementsByTagName('reference_2')[0]
                        refs.append(ref2.childNodes[0].data)
                        
                        if len(ref.getElementsByTagName('reference_3')) != 0:
                            ref3 = ref.getElementsByTagName('reference_3')[0]
                            refs.append(ref3.childNodes[0].data)
                            
                            if len(ref.getElementsByTagName('reference_4')) != 0:
                                ref4 = ref.getElementsByTagName('reference_4')[0]
                                refs.append(ref4.childNodes[0].data)
                                
                                if len(ref.getElementsByTagName('reference_5')) != 0:
                                    ref5 = ref.getElementsByTagName('reference_5')[0]
                                    refs.append(ref5.childNodes[0].data)
                reflist.append(refs)
            else:
                reflist.append('N/A')
            
            #HF298 list
            if hf298s.length != 0: # if <hf298> is not empty
                for hf in hf298s:
                    hfvalues = []
                    hf298_1 = hf.getElementsByTagName('hf298_1')[0] # first, take <hf298_1>
                    hfvalue_1 = hf298_1.childNodes[0].data
                    hfvalues.append(hfvalue_1)
                    
                    if len(hf.getElementsByTagName('hf298_2')) != 0: # if there is <hf298_2>
                        hf298_2 = hf.getElementsByTagName('hf298_2')[0] # add <hf298_2>
                        hfvalue_2 = hf298_2.childNodes[0].data
                        hfvalues.append(hfvalue_2)
                        
                        if len(hf.getElementsByTagName('hf298_3')) != 0: # if there is <hf298_3>
                            hf298_3 = hf.getElementsByTagName('hf298_3')[0] # add <hf298_3>
                            hfvalue_3 = hf298_3.childNodes[0].data 
                            hfvalues.append(hfvalue_3)
                            
                            if len(hf.getElementsByTagName('hf298_4')) != 0: # if there is <hf298_4>
                                hf298_4 = hf.getElementsByTagName('hf298_4')[0] # add <hf298_4>
                                hfvalue_4 = hf298_4.childNodes[0].data
                                hfvalues.append(hfvalue_4)
                                
                                if len(hf.getElementsByTagName('hf298_5')) != 0: # if there is <hf298_5>
                                    hf298_5 = hf.getElementsByTagName('hf298_5')[0] # add <hf298_5>
                                    hfvalue_5 = hf298_5.childNodes[0].data
                                    hfvalues.append(hfvalue_5)
                    converted_units = ConvertUnits(hfvalues)
                    hfvalue_main = ExtractMinMain(converted_units)[0]
                    min_pos = ExtractMinMain(converted_units)[1]
                    hfvalue = converted_units[min_pos] # this is the final HF value for <hf298>s, but <additional_information> also needs to be considered
                
            else:
                hfvalue = 'NaN'

            # check for HF298 value and reference in <additional_information>
            if additionals.length != 0: # additional information provided
                for adt in additionals:
                    hfvalues_adt1 = [] # collect all the hfs in <additional_information_1>
                    reflist_adt1 = [] # collect all the REFs in <additional_information_1>
                    adt1 = adt.getElementsByTagName('additional_information_1')[0]
                    adt_inf1 = adt1.childNodes[0].data
                    adt_inf1 = adt_inf1.upper()
                    count_hf2981 = adt_inf1.count('HF298=')
                    if count_hf2981 > 0: # HF298 provided in <additional_information_1>
                        hf_start_positions1 = FindSubStr('HF298=',adt_inf1,count_hf2981)
                        
                        hf_end_positions1 = []
                        ref_start_positions1 = []
                        for pos in range(len(hf_start_positions1)):
                            nearest_hf0 = adt_inf1.find('HF0=',hf_start_positions1[pos])
                            nearest_ref = adt_inf1.find('REF=',hf_start_positions1[pos])
                            if nearest_hf0 != -1 and nearest_ref != -1:
                                hf_end_positions1.append(min(nearest_hf0,nearest_ref))
                                ref_start_positions1.append(nearest_ref)
                            elif nearest_hf0 == -1 and nearest_ref != -1:# no HF0 after this hf298
                                hf_end_positions1.append(nearest_ref)
                                ref_start_positions1.append(nearest_ref)
                            elif nearest_hf0 == -1 and nearest_ref == -1:
                                hf_end_positions1.append(len(adt_inf1))
                                ref_start_positions1.append(-1)
                        for i in range(count_hf2981):
                            hfvalues_adt1.append(adt_inf1[hf_start_positions1[i]+6:hf_end_positions1[i]])
                            if i < count_hf2981-1:
                                reflist_adt1.append(adt_inf1[ref_start_positions1[i]:hf_start_positions1[i+1]])
                            else:
                                reflist_adt1.append(adt_inf1[ref_start_positions1[i]:])
                                
                        converted_units_adt1 = ConvertUnits(hfvalues_adt1)
                        hfvalue_main_adt1 = ExtractMinMain(converted_units_adt1)[0]
                        min_pos_adt1 = ExtractMinMain(converted_units_adt1)[1]
                        hfvalue_a1 = converted_units_adt1[min_pos_adt1] # this is the final HF value in adt_inf_1
                        
                        specie_ref_adt.append(reflist_adt1) # this is the list of refs in adt_inf_1
                        hfvalue_main_adt = hfvalue_main_adt1
                        hfvalue_adt = hfvalue_a1
                        
                    # look at the second adt_inf
                    if len(adt.getElementsByTagName('additional_information_2')) != 0:
                        hfvalues_adt2 = [] # collect all the hf values in <additional_information_2>
                        reflist_adt2 = [] # collect all the REFs in <additional_information_2>
                        adt2 = adt.getElementsByTagName('additional_information_2')[0]
                        adt_inf2 = adt2.childNodes[0].data
                        adt_inf2 = adt_inf2.upper()
                        count_hf2982 = adt_inf2.count('HF298=')
                        if count_hf2982 > 0: # there is HF298 in <additional_information_2>
                            hf_start_positions2 = FindSubStr('HF298=',adt_inf2,count_hf2982)
                            hf_end_positions2 = []
                            ref_start_positions2 = []
                            for pos in range(len(hf_start_positions2)):
                                nearest_hf0 = adt_inf2.find('HF0=',hf_start_positions2[pos])
                                nearest_ref = adt_inf2.find('REF=',hf_start_positions2[pos])
                                if nearest_hf0 != -1 and nearest_ref != -1:
                                    hf_end_positions2.append(min(nearest_hf0,nearest_ref))
                                    ref_start_positions2.append(nearest_ref)
                                elif nearest_hf0 == -1 and nearest_ref != -1:# no HF0 after this hf298
                                    hf_end_positions2.append(nearest_ref)
                                    ref_start_positions2.append(nearest_ref)
                                elif nearest_hf0 == -1 and nearest_ref == -1:
                                    hf_end_positions2.append(len(adt_inf1)-1)
                                    ref_start_positions1.append(-1)
                            for i in range(count_hf2982):
                                hfvalues_adt2.append(adt_inf2[hf_start_positions2[i]+6:hf_end_positions2[i]])
                                if i < count_hf2982-1:
                                    reflist_adt2.append(adt_inf2[ref_start_positions2[i]:hf_start_positions2[i+1]])
                                else:
                                    reflist_adt2.append(adt_inf2[ref_start_positions2[i]:])
                                    
                            converted_units_adt2 = ConvertUnits(hfvalues_adt2)
                            hfvalue_main_adt2 = ExtractMinMain(converted_units_adt2)[0]
                            min_pos_adt2 = ExtractMinMain(converted_units_adt2)[1]
                            hfvalue_a2 = converted_units_adt2[min_pos_adt2] # this is the final HF value for <hf298>s in adt_inf_2
                            specie_ref_adt.append(reflist_adt2) # this is the list of refs in adt_inf_2
                            
                            if hfvalue_main_adt2<hfvalue_main_adt1:
                                hfvalue_main_adt = hfvalue_main_adt2
                                hfvalue_adt = hfvalue_a2
    
                if hfvalue == 'NaN': #if <hf298> not provided, take the hf298 in <additional_inf>
                    hflist.append(hfvalue_adt)
                else: #If hf298 provided in both <hf298> and <adt_inf>, want the lowest one. If equal, take the one with uncertainty
                    if hfvalue_main == hfvalue_main_adt:
                        if ('+/-' in hfvalue) and ('+/-' not in hfvalue_adt):
                            hflist.append(hfvalue)
                        else:
                            hflist.append(hfvalue_adt)
                    elif hfvalue_main < hfvalue_main_adt:# if not equal, take the lower one 
                        hflist.append(hfvalue)
                    else:
                        hflist.append(hfvalue_adt)
            
            else: # additional information not provided
                hflist.append(hfvalue)
                specie_ref_adt.append(['No Additional Information Provided'])
                
            reflist_adt.append(specie_ref_adt)    
            
            for h in high_coef:
                specie_high.append(h.childNodes[0].data)
            highs.append(specie_high)    
            for l in low_coef:
                specie_low.append(l.childNodes[0].data)
            lows.append(specie_low)


print('There are ' + str(len(substance_l)) + ' gas phase species in the original dataset.')
print('Next, rewrite the substance names in a more standard way.')
print('=========================================================')



substancelist=[]
for thissub in substance_l:
    # substitute names by formulas
    thissub = re.sub('PHENYL-CN','C7H5N',thissub)
    thissub = re.sub('4-ETHENYL - PHENYL-1-VINYL','C6H4(C2H)C2H2',thissub)
    thissub = re.sub('CYCLOBUTADIENE','C4H4',thissub) 
    thissub = re.sub('BENZALDEHYDE','C7H6O',thissub)
    thissub = re.sub('S-OH','HOS',thissub)
    thissub = re.sub('N2D2-CIS','N2D2',thissub)
               
    # check cases 
    thissub = re.sub('AL','Al',thissub)
    thissub = re.sub('CL','Cl',thissub)
    thissub = re.sub('CA','Ca',thissub)
    thissub = re.sub('AR','Ar',thissub)
    thissub = re.sub('BR','Br',thissub)
    thissub = re.sub('BA','Ba',thissub)
    thissub = re.sub('FE','Fe',thissub)
    thissub = re.sub('HG','Hg',thissub)
    thissub = re.sub('MG','Mg',thissub)
    thissub = re.sub('GE','Ge',thissub)
    thissub = re.sub('MO','Mo',thissub)
    thissub = re.sub('CR','Cr',thissub)
    thissub = re.sub('HE','He',thissub)
    thissub = re.sub('NE','Ne',thissub)
    thissub = re.sub('XE','Xe',thissub)
    thissub = re.sub('ZN','Zn',thissub)
    thissub = re.sub('CDO','CdO',thissub)  
    # SN can be two elements (S & N) or Sn
    thissub = re.sub('SNCl4','SnCl4',thissub)  
    thissub = re.sub('SNH3','SnH3',thissub)
    thissub = re.sub('SNH4','SnH4',thissub)
    
    substancelist.append(thissub)



# clear form: substance formula only, phase not needed
print('Clear the Substance Formulas:')
substances = []
for i in range(len(substancelist)):
    g_pos = substancelist[i].find('(G') 
    if g_pos != -1:
        print(i, substancelist[i])
        substances.append(substancelist[i][:g_pos])
    else:
        substances.append(substancelist[i])


# clear form: symbols not needed
for i in range(len(substancelist)):
    h_pos = substancelist[i].find('##!!##') 
    if h_pos != -1:
        print(i, substancelist[i])
        findspace1 = substances[i].find(' ')
        findspace2 = substances[i].find(' ',7)
        substances[i] = substances[i][findspace1+1:findspace2]


# clear form: symbols not needed
for i in range(len(substancelist)):
    if substancelist[i][0] == '?':
        print(i, substancelist[i])
        findspace1 = substances[i].find(' ')
        findspace2 = substances[i].find(' ',4)
        substances[i] = substances[i][findspace1+1:findspace2]


# now the first word in every 'substance' is the rough formula
for i in range(len(substancelist)):
    space_pos = substances[i].find(' ')
    if space_pos != -1:
        substances[i] = substances[i][:space_pos]


for i in range(len(substancelist)):
    star_pos = substances[i].find('*')
    if star_pos != -1:
        print(i, substances[i])
        new_sub = ''
        for j in range(len(substances[i])):
            if j != star_pos:
                new_sub = new_sub + substances[i][j]
        substances[i] = new_sub
print('=========================================================')


for i in range(len(substancelist)):
    if '?' in substances[i] or '!' in substances[i] or '#' in substances[i] or '*' in substances[i]:
        print('ERROR: there is still error with MAGNESIUM, position: ' + str(i))
print('Number signs, exclamation marks, question marks, and asterisks are cleared! \n=========================================================')


# find errors in formulas of substances consist of MAGNESIUM
for i in range(len(substancelist)):
    mg_pos = substances[i].find('MAGNeSIUM')
    if mg_pos != -1:
        print('MAGNESIUM formula errors: ' + str(i) + '  ' + substancelist[i])
print('=========================================================')

substances[1032] = 'MgCl+'
substances[1033] = 'MgClF'
substances[1036] = 'MgF+'
substances[1042] = 'MgN'
substances[1048] = 'Mg2'

for i in range(len(substancelist)):
    if 'MAGNeSIUM' in substances[i]:
        print('ERROR: there is still error with MAGNESIUM, position: ' + str(i))
print('MAGNESIUM formula errors corrected! \n=========================================================')


n1 = 0
print('ERRORS with the Isomers (find hyphens):')
# there are three types of errors with hyphen
for i in range(len(substancelist)):
    if '-' in substances[i]:
        if substances[i][1] == '-' and substances[i][-1] != '-':
            print(i,substances[i])
            substances[i] = substances[i][2:]
            n1 = n1 + 1
        elif substances[i][-2] == '-':
            print(i,substances[i])
            substances[i] = substances[i][:-2]
            n1 = n1 + 1


n2 = 0
for i in range(len(substancelist)):
    if ',' in substances[i]:
        print(i,substances[i])
        hyphen_pos = substances[i].find('-')
        substances[i] = substances[i][hyphen_pos+1:]
        n2 = n2 + 1
print('=========================================================')


print('ERRORS with the Chemical Bonds (find hyphens):')
n3 = 0
for i in range(len(substancelist)):
    if '-' in substances[i] and substances[i][-1] != '-':
        print(i,substances[i])
        hyphen_pos = substances[i].find('-')
        new_sub = ''
        for j in range(len(substances[i])):
            if j != hyphen_pos:
                new_sub = new_sub + substances[i][j]
        substances[i] = new_sub
        n3 = n3 + 1
print('=========================================================')


for i in range(len(substancelist)):
    if '-' in substances[i] and substances[i][-1] != '-':
        print('ERROR: there is still error with hyphens, position: ' + str(i))
n = n1 + n2 + n3
print(str(n) + ' hyphen errors corrected! \n=========================================================')


print('ERRORS with the Chemical Bonds (find equal signs):')
for i in range(len(substancelist)):
    if '=' in substances[i]:
        print(i,substances[i])
        hyphen_pos = substances[i].find('=')
        new_sub = ''
        for j in range(len(substances[i])):
            if j != hyphen_pos:
                new_sub = new_sub + substances[i][j]
        substances[i] = new_sub
print('=========================================================')



for i in range(len(substancelist)):
    if '=' in substances[i]:
        print('ERROR: there is still error with equal sign, position: ' + str(i))
n = n1 + n2 + n3
print('Equal sign errors corrected! \n=========================================================')


n = 0
deuterium_mixture_pos = []
for i in range(len(substancelist)):
    if 'D' in substances[i] or substances[i] == 'AIR' or substances[i] == 'JET' or substances[i] == 'ELECTRON':
        deuterium_mixture_pos.append(i)
        n = n + 1


all_subs = [substances[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

print('Now the number of substances is ' + str(len(all_subs)) + '.')
print('=========================================================')



# clear () in the formulas
print('following brackets need to be cleared:')
for i in range(len(all_subs)):
    if '(' in all_subs[i]:
        print(i,all_subs[i])
        b_pos = all_subs[i].find('(')
        if all_subs[i][b_pos+1].isdigit():
            all_subs[i] = all_subs[i][:b_pos]

for i in range(len(all_subs)):
    if '(' in all_subs[i]:
        all_subs[i] = Deconstruct(all_subs[i])


for i in range(len(all_subs)):
    all_subs[i] = SumupElements(Addones(all_subs[i]))


# build rough result

for i in range(len(hflist)):
    hflist[i] = re.sub('NAN','NaN',hflist[i])
    if ((hflist[i] != 'NaN') and ('KJ' not in hflist[i])):
        hflist[i] = hflist[i] + ' KJ'


uncertainty = []
for i in range(len(hflist)):
    unctt = hflist[i].find('+')
    
    if unctt != -1 and hflist[i] != 'NaN':
        findkj = hflist[i].find('KJ')
        uncertainty.append(hflist[i][unctt+3:findkj])
        hflist[i] = hflist[i][:unctt]
        
    elif unctt == -1 and hflist[i] != 'NaN':
        findkj = hflist[i].find('KJ')
        hflist[i] = hflist[i][:findkj]
        uncertainty.append('NaN')
        
    else:
        uncertainty.append('NaN')
        

all_hf298 = [hflist[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_uncertainty = [uncertainty[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_high = [highs[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_low = [lows[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_ref_codes = [ref_codes[i] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_ref = [str(reflist[i])[2:-2] for i in range(len(substances)) if i not in deuterium_mixture_pos]

all_ref_adt = [str(reflist_adt[i])[1:-1] for i in range(len(substances)) if i not in deuterium_mixture_pos]


for i in range(len(all_ref_adt)):
    if all_ref_adt[i] == "['']":
        all_ref_adt[i] = "['No Additional Information Provided']"

f_data = open("rough_database.csv", "w", encoding = 'utf-8').close()
f_data = open("rough_database.csv", "a", encoding = 'utf-8')
f_data.write('Chemical_Substance,HF298 (kJ/mol),Uncertainty (kJ/mol),Coefficient_1 (High-temperature; T in Kelvin),Coefficient_2 (High-temperature; T in Kelvin),Coefficient_3 (High-temperature; T in Kelvin),Coefficient_4 (High-temperature; T in Kelvin),Coefficient_5 (High-temperature; T in Kelvin),Coefficient_6 (High-temperature; T in Kelvin),Coefficient_7 (High-temperature; T in Kelvin),Coefficient_1 (Low-temperature; T in Kelvin),Coefficient_2 (Low-temperature; T in Kelvin),Coefficient_3 (Low-temperature; T in Kelvin),Coefficient_4 (Low-temperature; T in Kelvin),Coefficient_5 (Low-temperature; T in Kelvin),Coefficient_6 (Low-temperature; T in Kelvin),Coefficient_7 (Low-temperature; T in Kelvin),Reference_Code,References,References_for_Additional_Information'+ '\n')


#using comma to seperate the grids, so need to check ','s in the reference texts
for i in range (len(all_ref_adt)):
    all_ref_adt[i] = re.sub(';','',all_ref_adt[i])
    all_ref_adt[i] = re.sub(',',';',all_ref_adt[i])
    all_ref_adt[i] = re.sub("'",'',all_ref_adt[i])


for i in range (len(all_ref)):
    all_ref[i] = re.sub(',',' ',all_ref[i])
    all_ref[i] = re.sub("'",'',all_ref[i])


for i in range(len(all_subs)): # air not needed
    f_data.write(all_subs[i] + ',' + all_hf298[i] + ',' + all_uncertainty[i] + ',' + all_high[i][0] + ',' + all_high[i][1] + ',' + all_high[i][2] + ',' + all_high[i][3] + ',' + all_high[i][4] + ',' + all_high[i][5] + ',' + all_high[i][6] + ',' + all_low[i][0] + ',' + all_low[i][1] + ',' + all_low[i][2] + ',' + all_low[i][3] + ',' + all_low[i][4] + ',' + all_low[i][5] + ',' + all_low[i][6] + ',' + all_ref_codes[i] + ',' + all_ref[i] + ',' + all_ref_adt[i][1:-1] + '\n')

f_data.close()
print('=========================================================\nThe rough result file rough_database.csv has been built.')



df = pd.read_csv('rough_database.csv',sep=',')

repeated_subs = df['Chemical_Substance']
duplicated = []


for i in range(1,len(repeated_subs)):
    if repeated_subs[i] == repeated_subs[i-1] and repeated_subs[i] not in duplicated:
        duplicated.append(repeated_subs[i])

deleted_rows = []
all_nan = []
for sub in duplicated:
    df_for_sub = df[df['Chemical_Substance'] == sub] # extract all rows for this substance
    rownames = df_for_sub.index.tolist()
    HFvalues_origin = df_for_sub['HF298 (kJ/mol)'].tolist()
    rows_to_compare = []
    hf_to_compare = []
    if all(i != i for i in HFvalues_origin): # if all the hfvalues of a substances are float 'nan', delete them all, and add them into a seperate list
        for i in rownames:
            deleted_rows.append(i)
            all_nan.append(i)
    else:
        for i in range(len(HFvalues_origin)):
            if str(HFvalues_origin[i]) != 'nan':# delete the lines with nan hf298. 
                rows_to_compare.append(rownames[i])
                hf_to_compare.append(HFvalues_origin[i])
        min_pos = rows_to_compare[hf_to_compare.index(min(hf_to_compare))] #this is the row that should be kept (lowest hf298)
        for i in rownames:
            if i != min_pos:
                deleted_rows.append(i)

df_new = df.drop(index = deleted_rows, inplace = False)

df_new.to_csv('Result_Database.csv', index = False, sep = ',', na_rep = 'NaN')
print('=========================================================\nThe result file Result_Database.csv has been built.')



df.iloc[all_nan].to_csv('Isomers_NaN_HF298.csv', index = False, sep = ',', na_rep = 'NaN')
csv_reader = csv.reader(open("Isomers_NaN_HF298.csv"))

nan_list=[]
for row in csv_reader:
    nan_list.append(row)

nan_list.remove(nan_list[0])

pos_nan = []
for i in nan_list:
    pos_nan.append(highs.index(i[3:10]))

original_names = []
for i in pos_nan:
    original_names.append(substance_l[i])

df_nan = df.iloc[all_nan]

df_nan.insert(1,'Original_Name_Structure',original_names)

df_nan.to_csv('Isomers_NaN_HF298.csv', index = False, sep = ',', na_rep = 'NaN')
print('=========================================================\nThe file for all the substances with no HF298 Isomers_NaN_HF298.csv has been built.')

# build ref code file
f_code = open("Ref_Code.csv", "w", encoding = 'utf-8').close()

f_code = open("Ref_Code.csv", "a", encoding = 'utf-8')
f_code.write('Reference codes,DATE origin' + '\n')
f_code.write('A,ARGONNE NAT.LABS. \nATcT A,Branko Ruscic; unpublished results from Active Thermochemical Tables v.1.25 using the Core (Argonne) Thermochemical Network v. 1.049 (May 2005). \nB,Ihsan Barin database. \n')
f_code.write('CODA,CODATA Tables. \nD,Delaware University. \nF,THERGAS calculations. \nIU,IUPAC data. \nJ,JANAF tables. \nG(L),NASA Glen(former Lewis) Research Center. \nP,Thermodynamic Research Center (Formerly American Petroleum Institute). \n')
f_code.write('R or Rus or TPIS,Russian Tables (TSIV/TPIS); Gurvich. \nS,Louisiana State University (LSU). \nT,Technion-Israel Inst. Technology. \nTT,New HF298 adjusted on old polynomial.')

print('=========================================================\nThe file for all the source codes Ref_Code.csv has been built.')
f_code.close()
