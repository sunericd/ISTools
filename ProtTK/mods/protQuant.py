from csv import reader
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression
import mods.saveOutput as saveOutput

def protQuant (standards_file, protein_file, N_mass=25, regression_style='Linear'):
    '''
    standards_file = csv or tsv file where first row are the standard concentrations and other rows are fluoresence
    protein_file = csv or tsv file where rows are fluoresences
    N_mass = FLOAT: amount of protein wanted (in same units as in input)
    regression_style = STRING: type of regression to be used: 'Linear', 
    '''
    try:
        float_test = float(N_mass)
    except:
        print ("Error: The 'N mass' needs to be an integer or decimal value!")
        return
    
    # Reading files
    standard_concs, standard_fluors = readStandards(standards_file)
    prot_fluors = readProteins(protein_file)
    
    # Averaging multiple samples...
    standard_fluors = averageRows(standard_fluors)
    prot_fluors = averageRows(prot_fluors)
    
    # Finding regression...
    prot_concs = regressor(standard_concs, standard_fluors, prot_fluors, regression_style)
    
    # Finding N volume needed for N mass...
    N_vols = []
    for prot_conc in prot_concs:
        try:
            N_vol = float(N_mass)/prot_conc
        except:
            print ("Error: Protein Concentration regresses to zero!")
            return
        N_vols.append(N_vol)
    
    out_text = str(prot_concs)
    for item in N_vols:
        out_text=out_text+("%s\n" % item)
        
    saveOutput.saveData(out_text, "ProtQuant")

def readStandards (filename):
    standards_concs = []
    standard_fluors = []
    
    # Reading standards...
    with open(filename) as file:
        if '.tsv' in filename:
            i = 0
            for line in reader(file, delimiter="\t"):
                if i == 0:
                    standard_concs = []
                    for item in line:
                        standard_concs.append(float(item))
                else:
                    standard_fluors.append(line)
                    # Raises error if number of cols doesn't match
                    if len(standard_fluors[i-1]) != len(standard_concs):
                        raise Exception('Please make sure that the input is rectangular.')
                i += 1
            
        elif '.csv' in filename:
            i = 0
            for line in reader(file, delimiter=","):
                if i == 0:
                    standard_concs = []
                    for item in line:
                        standard_concs.append(float(item))
                else:
                    standard_fluors.append(line)
                    if len(standard_fluors[i-1]) != len(standard_concs):
                        raise Exception('Please make the standard file rectangular.')
                i += 1
        else:
            raise Exception('Please input a tsv or csv file for Standards!')
    
    return (standard_concs, standard_fluors)

def readProteins (filename):
    prot_fluors = []
    
    # Reading standards...
    with open(filename) as file:
        if '.tsv' in filename:
            for line in reader(file, delimiter="\t"):
                prot_fluors.append(line)
        elif '.csv' in filename:
            for line in reader(file, delimiter=","):
                prot_fluors.append(line)
        else:
            raise Exception('Please input a tsv or csv file for Proteins!')
            
    if not all(len(prot_fluor) == len(prot_fluors[0]) for prot_fluor in prot_fluors):
        raise Exception('Please make the protein file rectangular.')
    
    return (prot_fluors)

def averageRows (list_of_lists):
    float_lists = []
    for x in list_of_lists:
        float_x = list(map(float, x))
        float_lists.append(float_x)
    avg = [float(sum(col))/len(col) for col in zip(*float_lists)]
    return (avg)

def regressor(standard_concs, standard_fluors, prot_fluors, regression_style):
    # Reshaping...
    standard_fluors = np.array(standard_fluors)
    standard_fluors = standard_fluors.reshape(len(standard_fluors), 1)
    prot_fluors = np.array(prot_fluors)
    prot_fluors = prot_fluors.reshape(len(prot_fluors), 1)
    # Finding the model...
    if regression_style == 'Linear':
        regr = LinearRegression()
    else:
        print ('Defaulting to Linear Regression...')
        regr = linear_model.LinearRegression()
    # Training
    regr.fit(standard_fluors, standard_concs)
    # Testing/Predicting
    prot_concs = regr.predict(prot_fluors)
    return (prot_concs)


#protQuant('test_standards.csv', 'test_proteins.csv', 25, 'Linear')