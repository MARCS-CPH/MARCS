import pandas as pd
import numpy as np
class MARCS_model():
    def __init__(self, file):
        self. lines = file.readlines()
        self.data = {}
    def general_info(self):
    #read general info and stores them in a dictionary
            i = 0
            
            lines = self.lines
            
            for line in lines:
                
                words = line.strip().split(' ')
                i+=1
                for j in range(0,len(words)):
                    #if words[j].isalpha():
                        
                        if words[j] == 'TEFF=':
                            
                            self.data['TEFF'] = float(words[j+1])

                        if words[j] == 'G=':
                            
                            self.data['LOGG'] = float(words[j+2])

                        if words[j] == 'SIGMA*TEFF**4)':
                            
                            self.data['Total Flux'] = float(words[j+2])

                        if words[j] == 'C/O':
                            
                            word = words[j+3]
                            
                            self.data['C/O'] = float(word)
                        if words[j] == 'GRAVITY':
                            x=(words[j+5])
                            
                            self.data['acceleration of gravity'] = float(x)
                        if len(words)>2:
                            if words[1] == 'NORMAL' or words[2]=='END':

                                print('file ended')
                            

            return self.data

    def thermo(self):
        #reads 'Correction from last iteration', 
        #'Model atmosphere' and ' Thermodynamical quantities and convectionCGS UNITS)'
        
        lines = self.lines
       
        constrain = 0
        df_list = list()
        obs_list = list()
        
        counter = 0
        # Strips the newline character
        
        for line in lines:
            my_string= line.strip().split(' ')
            if my_string[0] == 'P':
                break
            if constrain == 0:
                if my_string[0] == 'K':
                    col_names = list()
                    
                    for elm in my_string:
                        if len(elm) >0 and elm !='GEOM.' and elm != 'VEL.':
                            
                            col_names.append(elm)
                    constrain = 1
                   
                else:
                    continue
                   
            else:
                if len(my_string) > 1 and my_string[0]!='K':
                    obs_row = list()
                    for row in my_string:
                        if len(row) >0 :
                            #print(row)
                            obs_row.append(float(row))
                    obs_list.append(obs_row)
                    
                else:
                    df_list.append(pd.DataFrame(obs_list, columns = col_names).iloc[: , :-1])
                    obs_list=list()
                    col_names=list()
                    constrain = 0
            counter +=1
        
        return df_list 
     
    def partial_pressures(self):
        lines = self.lines
        constrain = 0
        df_list = list()
        obs_list = list()
        
        counter = 0
        i=0
        for line in lines:
            my_string= line.strip().split(' ')
            if constrain == 0:
                if my_string[0] == 'P':
                    col_names = list()
                    next_line = lines[i+1]
                    next_string = next_line.strip().split(' ')
                    if next_string[0] =='22':
                        break
                    else:
                        for elm in next_string:
                            if elm != '':
                                col_names.append(elm)
                        print(col_names)
                        constrain = 1
                        i +=1
                   
                else:
                    i +=1
                    continue
            else:
                if len(my_string) > 1 and my_string[0]!='P':
                    obs_row = list()
                    next_line = lines[i+1]
                    next_string = next_line.strip().split(' ')
                    for row in next_string:
                        if len(row) >0 :
                            obs_row.append(float(row))
                    #print(obs_row)
                    if len(row) >0:    
                        del obs_row[0]
                        obs_list.append(obs_row)
                    i+=1
                    
                else:
                    df_list.append(pd.DataFrame(obs_list, columns = col_names))
                    obs_list=list()
                    col_names=list()
                    constrain = 0
                    i+=1
            
            
        return df_list

input_file = open('converge.in', 'r')
lines = input_file.readlines()
filenames = []
for line in lines:
    filenames.append(line.strip().split('\n'))

old = filenames[0][0]
new = filenames[1][0]

df_old = MARCS_model(open(old, 'r')).thermo()
df_new = MARCS_model(open(new, 'r')).thermo()

delta_T = abs(df_old[1]['T']-df_new[1]['T'])/df_old[1]['T']
#delta_T = abs(df_old[1]['T']-df_new[1]['T']) / df_old[1]['T']
print(delta_T.max())
with open('check_conv.txt', 'w') as f:
    f.write(str(delta_T.max()))