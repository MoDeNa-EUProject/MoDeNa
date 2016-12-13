
import sys
import time

# The following lists are required to replace built-in Predici functions by corresponding fluent functions

predici_command_list=['gettemp','getmass','getvol','getco','cpol','getmmlow']
predici_command_arg_list=[0,0,0,2,2,2]
predici_command_type_list=['double','double','double','double','double','double']
predici_command_replace_list=['my_temp','my_cell_mass','my_cell_volume','my_conc[MY_NEW_INDEX]','my_conc[MY_NEW_INDEX]','my_molma[MY_NEW_INDEX]']
fluent_command_replace_list=['C_T(c,t)','C_VOLUME(c,t)*C_R(c,t)','C_VOLUME(c,t)','my_conc[MY_NEW_INDEX]','C_UDSI(c,t,MY_NEW_INDEX)','my_molma[MY_NEW_INDEX]']
fluent_command_array_size_list=['0','0','0','40','40','40']
target_list=['Openfoam','Fluent']

# Function types

func_type_list=['void','bool','int','long','double','float','char']
 
#===================================================================================================================================================================
# Definition of a method through a class which applied to a string and which has a substring list as argument.
# The method returns the substring in the list which is found first together with the position of the first occurence of the substring in the string 
        
def firststring(stringname,string_list):
    myindex = -1
    myname = ''
    for index in range(len(string_list)):
        index_sstr = stringname.find(string_list[index])
        if index_sstr >= 0:
            if myindex == -1:
                myindex = index_sstr
                myname = string_list[index]
            else:
                if myindex > index_sstr:
                    myindex = index_sstr
                    myname = string_list[index]
    return myname, myindex      
                    
#================================================================================================================================================================================
#
# Routine tests if a substring can be a function name
# Therefore it checks wether the characters before and behind substring in string are suitable
#
# Arguments:
# str1: String in which substring is checked as a possible function name
# substr1: The substring in string which is checked as possible function name
# i1: The position of the first character of substring in string (this is required since substr1 could occur several times


def check_func_name(str1,sub_str1,i1):
    stra = str1[:i1]
    strb = str1[i1+len(sub_str1):]

    stra = stra.replace('\n',' ') 
    stra = stra.replace('\t',' ') 
    strb = strb.replace('\t',' ')
    strb = strb.replace('\n',' ')
    strb = strb.lstrip() 

    if len(stra) == 0:
        icheck1 = 1
    elif stra[len(stra)-1:len(stra)] == ' ':
        icheck1 = 1
    else:
        icheck1 = 0

    if len(strb) == 0:
        icheck2 = 1
    elif strb[0:1] == '(':
        icheck2 = 1
    else:
        icheck2 = 0
                        
    if (icheck1 == 1) and (icheck2 == 1):
        return 1
    else:
        return 0
    return 0


#================================================================================================================================================================================
# Routine finds a corresponding bracket in a string
# 
# Arguments:
# str1: String in which correspoinding bracket is searched for
# str2: Starting bracket type
# i1:   Index of the starting bracket in the string

def corresponding_bracket_index(str1,str2,i1):
    index = i1+1
    if str2 not in ['(','{','[']:
        print 'Error in function corresponding_bracket_index. str2 is not an opening bracket'
        raise SystemExit
    else:
        if str2 == '(':
            str3 = ')'
        if str2 == '{':
            str3 = '}'
        if str2 == '[':
            str3 = ']'

    count_bracket = 1
    while index < len(str1) and count_bracket > 0:
       while str1[index:].find(str2)>=0 or str1[index:].find(str3)>=0:
           index_b1 = str1[index:].find(str2)
           index_b2 = str1[index:].find(str3)
           if (index_b1 == -1) and (index_b2 >=0):
               count_bracket -= 1
               index = index+index_b2+1
               if count_bracket <= 0:
                   return index-1
               else:
                   print 'Unbalanced bracket'
                   raise SystemExit
           if (index_b2 == -1) and (index_b1 >=0):
               print 'Unbalanced bracket'
               raise SystemExit 
           if (index_b1 < index_b2):
               count_bracket += 1
               index = index+index_b1+1
           elif (index_b2 < index_b1):
               count_bracket -= 1
               index = index+index_b2+1
               if count_bracket <= 0:
                   return index-1
    print 'Error in function corresponding_bracket_index. Wrong syntax of interpreted source code'
    raise SystemExit    
    return index-1

#==================================================================================================================================
# Function converts a function of type func_type into a void function. It removes therefore all return values 

def func_2_voidfunc(code_string,func_type):
    itype = len(func_type)
    while code_string.find(func_type) >= 0:
        ibool = code_string.find(func_type)
        temp_string1 = code_string[:ibool+itype]
        temp_string1 = temp_string1.replace(func_type,'void')
        temp_string2 = code_string[ibool+itype:]
        first_string, first_index = firststring(temp_string2,['(',';','='])
        if first_string == '(':
            # Function definition - Find enclosing curlibrackets of function body and remove return variable after return statement
            temp_string1 += temp_string2[:first_index]
            temp_string2 = temp_string2[first_index:]
            sc_index = temp_string2.find('{')
            ec_index = corresponding_bracket_index(temp_string2,'{',sc_index)
            temp_string1 += temp_string2[:sc_index]
            temp_string3 = temp_string2[ec_index:]
            temp_string2 = temp_string2[sc_index:ec_index]
            while temp_string2.find('return') >= 0:
                index1 = temp_string2.find('return')
                temp_string1 += temp_string2[:index1+6]
                temp_string2 = temp_string2[index1+6:]
                if temp_string2[0] == ' ':
                    index2 = temp_string2.find(';')
                    temp_string1 += temp_string2[index2]
                    temp_string2 = temp_string2[index2+1:]
            temp_string1 += temp_string2+temp_string3
            temp_string2 = ''
            temp_string3 = ''
        code_string = temp_string1+temp_string2

    return code_string

#================================================================================================================================================================================
# Routine identifies whether a function occurs in a string or not
# It returns the position of the first character of the function name in the string and the position of the closing bracket ')'. The closing bracket must be in the same string!
# If the function is not found in the string it returns -1,-1
#
# Arguments: 
# str1: String in which the function name is searched
# str2: Function name
# 
# Results:
# index_start: Position of the first character of the function name in the string
# index_end: Position of the closing bracket ')' in the string

def my_find_func( str1, str2 ):
    index4 = str1.find(str2)

    if index4  == -1:
        return -1,-1
 
    bindex = str1.find('(',index4+len(str2)) 
    if bindex == -1:
        return -1,-1
    else:
        str3 = str1[index4+len(str2):bindex]
        str3 = str3.lstrip()
        if len(str3) > 1:
            return -1,-1  

    pchar=''
    if index4 > 0: 
        pchar = str1[index4-1]     
        if pchar not in [' ','=','>','<','*','+','-','/','(',',']:
            return -1,-1
    index_start = index4
    index4 = bindex
    index5 = str1.find(')',index4,)
    while str1.find('(',index4+1) > -1 and str1.find('(',index4+1) < index5:
        index5 = str1.find(')',index5+1)
        index4_ = str1.find('(',index4+1)
    index_end = index5
   
    return index_start, index_end



    
#================================================================================================================================================================================
# This functions returns the argument list of a function call as a string
#
# Arguments:
# str: The string which contains the function call together with the argument list
# i1: The position of the first character of the function name in the string str
# i2: The position of the closing bracket ')' of the argument list
#
# Results:
# arg_string: string of the argument list of a function including ',' separators
#
# Basically the function checks for the opening bracket '(' and returns the substring between the position of '(' and the corresponding closing bracket ')'


def get_arg_list(str,i1,i2):
    indexb = str.find('(',i1,i2)

    arg_string = str[indexb:i2]

    return arg_string




#================================================================================================================================================================================
# This function splits a string 'str1' into a list of strings. The split string is 'str2'. From the created list the entry 'i1' is returned
#
# Arguments:
# str1: The string which is split
# str2: The separator string
# i1: The entry of the generated list that shall be returned
#
# Results:
# arg_string: The string of the i1 argument in the argument list
    
def get_list_entry(str1,str2,i1):
    arg_list = str1.split(str2)
    arg_string = arg_list[i1-1]
    
    return arg_string

#================================================================================================================================================================================  
# This function deletes the body of a C or C++ function (definition)
#
# Arguments:
# my_list: Contains the lines of the C/C++ code
# i1: Is the line in which the function starts with the function name
#
# Results:
# None
#
# Everything furtheron is deleted until a ';' or a '{' sign are found. Then there are two scenarios
# (1) If a ';' is found first it is also deleted and the deletion then stops.
# (2) If a '{' is found first everything until and including the corresponding '}' is deleted. Then the deletion process stops

def delete_func(code_string):
    my_delete_list =['exp','power']

    new_code = ''
    while len(code_string) > 0:
        functype,itype = firststring(code_string,func_type_list)
        funcname,ifunc = firststring(code_string,my_delete_list)
        if itype < 0 or ifunc < 0:
            new_code += code_string
            code_string = ''
        else:
            if ifunc < itype:
                new_code += code_string[:itype]
                code_string = code_string[itype:]
            else:
                test_string = code_string[itype+len(functype):ifunc]
                test_string = test_string.strip()
                if test_string == '':
                    istart,iend = my_find_func(code_string,funcname)
                    new_code += code_string[:itype]
                    code_string = code_string[iend:]
                    mystring,myindex=firststring(code_string,[';','{'])
                    if mystring == ';':
                        code_string = code_string[myindex:]
                    else:
                       iend = corresponding_bracket_index(code_string,'{',myindex)
                       code_string = code_string[iend+1:]
                else:
                    new_code += code_string[:itype+len(functype)]
                    code_string = code_string[itype+len(functype):]
    return new_code

#=================================================================================
# Function finds the position of a the first for-loop in a string
#
# Arguments: 1
#
# str1 = The string in which the for loop is searched 
#
# Results: 1
#
# index = The position of the first for-loop in the string 
#

def find_loop(str1):
    index = 0
    while index < len(str1):
        index = str1[index:].find('for')+index
        if index >=0:
            str2 = str1[index+3:]
            str2 = str2.lstrip()
            if str2[0:1] == '(':
                if index == 0:
                    return index
                elif str1[index-1:index] in [' ',';','}','{']:
                    return index
                else:
                    index += 3
            else:
                index += 3
        else:
            return -1
        
    return -1

#================================================================================================
# Function modifies a for-loop in a string which is located at a given position in the string if 
# it contains an internal declaration of the loop counter:
# First it is checked whether an internal declaration exists in this for-loop.
# If yes, the internal declaration is removed and an external declaration is added in front of the loop
#
# Arguments: 2
# str1 = String which contains the for-loop
# i1 = Starting position of the for-loop in the string
#
# Results: 1
# str3 = String that is returned. If no internal counter declaration was found, str3 is the same as the input string str1

def remove_loop_int_declaration(str1,i1):
    index=i1+str1[i1:].find('(')
    str2 = str1[index+1:]
    str2 = str2.lstrip()
    if str2.find('int ') == 0:
        index_start = str1[index:].find('int ')+index
        index_end = index_start+4
        str2=str2[4:].lstrip()
        str2=str2[:str2.find('=')]
        str2=str2.rstrip()
        str3 = str1[:i1]+'int '+str2+';'+str1[i1:index_start]+str1[index_end:]
    else:
        str3 = str1

    return str3

#================================================================================================
# Function modifies all those for-loops in a string that contain an internal declaration of the 
# loop counter
#
# Arguments: 1
# str1 = String in which the for-loops shall be changed
#
# Results: 1
# new_code = Return string in which all for-loops which internal counter declaration have been c
#            changed

def remove_all_loop_int_dec(str1):
    new_code =''
    while len(str1) > 0:
        iloop = find_loop(str1)
        str2 = remove_loop_int_declaration(str1,iloop)
        new_code += str2[:iloop+3]
        str1 = str2[iloop+3:]
    return new_code         
        
#==================================================================================================
# function removes unnecessary semicolons
#
# Arguments: 1
# str1 = Input string which potentially contains unnecessary semicolons
#
# Results: 1
# str2 = Return string which contains only necessary semicolons

def remove_semicolon(str1):
    my_list_f1 = ['{','}']
    my_list_f2 = ['{']
    my_list_b1 = [';']
    my_list_b2 = ['}']
    str2 =''
    while len(str1) > 0:
        index = str1.find(';')
        if index < 0:
            str2 += str1
            str1 = ''
        else:
            temp_string = str1[index+1:].replace(' ','') 
            if temp_string[0] in my_list_b1:
                str3, delindex = firststring(str1[index+1:],my_list_b1)
                str1 = str1[:index+1]+str1[index+1+delindex+1:]
            else:
                temp_string=str1[:index+1].replace(' ','') 
                if temp_string[len(temp_string)-2] in my_list_f1:
                    if temp_string[len(temp_string)-2] not in my_list_f2:
                        temp_string = str1[index+1:].replace(' ','')
                        if temp_string[0] in my_list_b2:
                          str1 = str1[:index]+str1[index+1:]
                        else:
                          str2 += str1[:index+1]
                          str1 = str1[index+1:]   
                    else:
                        str1 = str1[:index]+str1[index+1:]
                else:
                    temp_string=str1[:index+1].replace(' ','')
                    i1 = len(temp_string)-5
                    i2 = len(temp_string)-1
                    temp_string=temp_string[i1:i2]
                    if  temp_string == 'else':
                       str1 = str1[:index]+str1[index+1:]
                    else:
                        str2 += str1[:index+1]
                        str1 = str1[index+1:]

    return str2

#==================================================================================================
# function changes pass by reference to c style constructs
#
# Arguments: 1
# str1 = Input string wich contains the code where passing by reference has to be replaced
#
# Results: 1
# str2 = Return string which contains the code where passing by reference has been replaced

def pass_by_reference_2_c(str1):
    index = 0
    searchstring='double& result1, double& result2'
    replacestring ='double* results'
    while index < len(str1):
        deltai = str1[index:].find(searchstring)
        if deltai < 0:
            index = len(str1)
        else:
            index = index+deltai
            icurli = index+str1[index:].find('{')
            isemi = index+str1[index:].find(';')
            if icurli >= 0 and (icurli < isemi or isemi < 0):
                sub0 = str1[index:]
                ibend = index+corresponding_bracket_index(sub0,'{',icurli-index)
                sub1 = str1[:icurli]
                sub2 = str1[icurli:ibend+1]
                sub3 = str1[ibend+1:]
                sub2 = sub2.replace('result1','results[0]')
                sub2 = sub2.replace('result2','results[1]')
                str1 = sub1+sub2+sub3
                index = len(sub1)+len(sub2)
            else:
                index +=len(searchstring)
    str2 =''
    str1 = str1.replace(searchstring,replacestring)
    while len(str1) > 0:
        index = str1.find('k1, k2, k1, k2')
        if index < 0:
            str2 += str1
            str1 = ''
        else:
            str2 += str1[:index-1]
            str1 = str1[index-1:]
            isemi = str1.find(';')
            str1 = 'k1, k2, results); k1 = results[0]; k2 = results[1];'+str1[isemi+1:]

    index = str2.find('void F')
    icurli = index + str2[index:].find('{')
    sub1 = str2[:icurli+1]
    sub2 = 'double results[2];'
    sub3 = str2[icurli+1:]
    str2 = sub1+sub2+sub3
    index = str2.find('void A')
    icurli = index + str2[index:].find('{')
    sub1 = str2[:icurli+1]
    sub2 = 'double results[2];'
    sub3 = str2[icurli+1:]
    str2 = sub1+sub2+sub3
    
    return str2

#=====================================================================================
# Reformating code

def reformat_code(lines):

    for index in range(len(lines)):
        new_string = lines[index]
        icom = new_string.find('//')
        if icom >= 0:
            lines[index]=new_string[0:icom]

    for index in range(len(lines)):
        new_string = lines[index]
        new_string = new_string.strip(' \t\r\n')
        new_string = new_string.replace('\t',' ')
        while new_string.find('  ') >= 0:
            new_string = new_string.replace('  ',' ')
        lines[index] = new_string

    code_string =''

    for index in range(len(lines)):
        code_string+=lines[index]
    return code_string


#===================================================================================
# Add new lines

def add_newlines(code_string):
    code_string = code_string.replace (';',';\n')
    code_string = code_string.replace ('{','\n{\n')
    code_string = code_string.replace ('}','}\n')
    
    return code_string

#===================================================================================
# Add F and A calling sequence

def add_call_FandA(mystring):
    mystring+='void FandA(double* x, double t, double* fall){int i=0; double fx[DIM]; double dfx[DIM]; F(x,t,fx); A(x,t,dfx); for (i=0; i< DIM; i++){fall[i]=fx[i];fall[i+DIM]=dfx[i];} return;}'
    return mystring
#==========================================================================================================================================================          
# Start of main program   

def create_args(kinetics_name):

    fin = open(kinetics_name, 'r')

    lines=fin.readlines()

    fin.close()

#====================================================================================
# Reformat code

    code_string = reformat_code(lines)

#====================================================================================
# Change bool functions to void function and remove return values
 
    code_string = func_2_voidfunc(code_string,'bool')

    code_string = delete_func(code_string)

#====================================================================================
# Find and modify int definitions in loops and transfer to C99 style
#

    code_string =  remove_all_loop_int_dec(code_string)

#====================================================================================
# Insert parameter initialization in void F function

    index = code_string.find('void F')

    index1 = index+code_string[index:].find('{')
    initcode ='if (my_global_init == -1) {SetGlobalIniParam();my_global_init=1;}'
    code_string = code_string[:index1+1]+initcode+code_string[index1+1:]

#=====================================================================================
# Remove unnecessary semicolons
#
    code_string = remove_semicolon(code_string)

#=====================================================================================
# Convert passing by reference to C style
#
    code_string = pass_by_reference_2_c(code_string)

    code_string = add_newlines(code_string)
  
    new_lines = code_string.split('\n')

#===================================================================================
# Include standard lib math.h in code in first line

    new_lines.insert(0,'#include <math.h>\n')
    new_lines.insert(1,'int my_global_init = -1;\n')

    str3 = ''
    for index in range(len(new_lines)): 
        str3+=new_lines[index]

# ---------------------- Create predici header file ----------------------------------

    predicicode = add_call_FandA(str3)

    predicicode = add_newlines(predicicode)

    # ---------------------- Create function code ----------------------------------------

    headerfilelist=['modena.h','math.h']

    codestring =''

    for index in range(len(headerfilelist)):
        codestring += '#include "'+headerfilelist[index]+'"\n'

    codestring += '''

void FandA(double* x, double t, double* fall);


void predici_kinetics
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    double time = inputs[0];
    double * inputParams = (double *) (&(inputs[1]));
    FandA(inputParams, time, outputs);
}
'''

    codestring += predicicode

    # ------------ Extract input list -------------------- 

    index = str3.find('XNAMES')
    strnames=str3[index:len(str3)]
    istart=strnames.find('{')
    iend=strnames.find('}')
    strnames=strnames[istart+1:iend]
    namelist = strnames.split(',')
    for index in range(len(namelist)):
        namelist[index] = namelist[index].replace('"',"'")

    namelist.insert(0,"'kineticTime'");
    returnDict = {}
    returnDict['Ccode'] = codestring
    returnDict['inputs'] = {}
    returnDict['outputs'] = {}
    returnDict['parameters'] = {}

    for index in range(len(namelist)):
        returnDict['inputs'][namelist[index]] = \
            {'min': -10, 'max': 9e99, 'argPos': index}
    namelist.pop(0);

    # ------------ Extract output list --------------------

    index = str3.find('XTYPES')
    strnames=str3[index:len(str3)]
    istart=strnames.find('{')
    iend=strnames.find('}')
    strnames=strnames[istart+1:iend]
    typelist = strnames.split(',')

    index_out = 0

    for index in range(len(namelist)):
        source_name='source_'+namelist[index].strip("'")
        returnDict['outputs'][source_name] = \
            {'min': -9e99, 'max': 9e99, 'argPos': index_out}
        index_out += 1

    for index in range(len(namelist)):
        dsource_name='dsource_'+namelist[index].strip("'")
        returnDict['outputs'][dsource_name] = \
            {'min': -9e99, 'max': 9e99, 'argPos': index_out}
        index_out += 1

    return returnDict

