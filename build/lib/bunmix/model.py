import ipywidgets as ipyw
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick
import bunmix.burr as br
import bunmix.input as bi

def get_Cfactor(X):
    ##### UNIT CONVERSION FROM EMU #####
    # Convert emu to emu/cm3
    iunits = X['iunits_widge'].value
    ounits = X['ounits_widge'].value
    volume = X['volume_widge'].value
    mass = X['mass_widge'].value

    Cfactor = 1

    if (iunits==0) & (ounits==1):
            Cfactor = 1/volume

    # Convert emu to emu/g
    if (iunits==0) & (ounits==2):
            Cfactor = 1/mass
    
    # Convert emu to Am2
    if (iunits==0) & (ounits==3):
            Cfactor = 0.001    
    
    # Convert emu to A/m
    if (iunits==0) & (ounits==4):
            Cfactor = 0.001/(volume/1E6)
            
    # Convert emu to Am2/kg
    if (iunits==0) & (ounits==5):
            Cfactor = 1/mass

    ##### UNIT CONVERSION FROM EMU/CM3 #####
    # Convert emu/cm3 to emu
    if (iunits==1) & (ounits==0):
            Cfactor = volume
    
    
    # Convert emu/cm3 to emu/g
    if (iunits==1) & (ounits==2):
            Cfactor = volume/mass
    
    # Convert emu/cm3 to Am2
    if (iunits==1) & (ounits==3):
            Cfactor = volume/0.001    
    
    # Convert emu/cm3 to A/m
    if (iunits==1) & (ounits==4):
            Cfactor = volume*0.001
            
    # Convert emu/cm3 to Am2/kg
    if (iunits==1) & (ounits==5):
            Cfactor = volume/mass

    ##### UNIT CONVERSION FROM EMU/G ###############################
    # Convert emu/g to emu
    if (iunits==2) & (ounits==0):
            Cfactor = mass
    
    # Convert emu/g to emu/cm3
    if (iunits==2) & (ounits==1):
            Cfactor = mass/volume
    
    # Convert emu/g to Am2
    if (iunits==2) & (ounits==3):
            Cfactor = mass*0.001    
    
    # Convert emu/g to A/m
    if (iunits==2) & (ounits==4):
            Cfactor = mass/volume*1000
            
    # Convert emu/g to Am2/kg
    if (iunits==2) & (ounits==5):
            error = 0
            Cfactor = 1

            
    ##### UNIT CONVERSION FROM Am2 ###############################
    # Convert Am2 to emu
    if (iunits==3) & (ounits==0):
            error = 0
            Cfactor = 1000
    
    # Convert Am2 to emu/cm3
    if (iunits==3) & (ounits==1):
            Cfactor = 1000/volume
    
    # Convert Am2 to emu/g
    if (iunits==3) & (ounits==2):
            Cfactor = 1000/mass    
    
    # Convert Am2 to A/m
    if (iunits==3) & (ounits==4):
            Cfactor = 1/(volume/1E6)
            
    # Convert Am2 to Am2/kg
    if (iunits==3) & (ounits==5):
            Cfactor = 1/(mass/1000)

    ##### UNIT CONVERSION FROM A/m ###############################
    # Convert A/m to emu
    if (iunits==4) & (ounits==0):
            Cfactor = 1000*volume
    
    # Convert A/m to emu/cm3
    if (iunits==4) & (ounits==1):
            error = 0
            Cfactor = 1000
    
    # Convert A/m to emu/g
    if (iunits==4) & (ounits==2):
            Cfactor = 1000*volume/mass    
    
    # Convert A/m to Am2
    if (iunits==4) & (ounits==3):
            Cfactor = volume/1E6
            
    # Convert A/m to Am2/kg
    if (iunits==4) & (ounits==5):
            Cfactor = (volume/1E6)/(mass/1000)

    ##### UNIT CONVERSION FROM Am2/kg ###############################
    # Convert Am2/kg to emu
    if (iunits==5) & (ounits==0):
            Cfactor = mass
    
    # Convert Am2/kg to emu/cm3
    if (iunits==5) & (ounits==1):
            Cfactor = mass/volume
    
    # Convert Am2/kg to emu/g
    if (iunits==5) & (ounits==2):
            Cfactor = 1    
    
    # Convert Am2/kg to Am2
    if (iunits==5) & (ounits==3):
            Cfactor = mass/1E3
            
    # Convert Am2/kg to A/m
    if (iunits==5) & (ounits==4):
            Cfactor = (mass/1000)/(volume/1E6)

    X['Cfactor'] = Cfactor

    return X

def get_params(X):
    
    X = get_Cfactor(X)

    nB_widge = ipyw.Dropdown(style={'description_width': 'initial'},
    options=[('1',0),('2',1),('3',2),('4', 3)],
    value=0,description='Number of components:')

    nTune_widge = ipyw.IntText(style={'description_width': 'initial'},
    value=2000,
    description='Number of tuning steps [>100]:',
    disabled=False)
    
    nSample_widge = ipyw.IntText(style={'description_width': 'initial'},
    value=2000,
    description='Number of sampling steps [>100]:',
    disabled=False)
    
    params_widge = ipyw.VBox([nB_widge,nTune_widge,nSample_widge])

    params_nest = ipyw.Tab()
    params_nest.children = [params_widge]
    params_nest.set_title(0, 'Model Parameters')
    display(params_nest)
    
    X['nB_widge'] = nB_widge
    X['nTune_widge'] = nTune_widge
    X['nSample_widge'] = nSample_widge
    
    return X

def get_model(X):
    
    #Number of components
    nB = X['nB_widge'].value+1
    
    #load data
    X = bi.load_data(X)
    x = X['B']
    y = X['M']

    if X['type_widge'].value==0:
        x1 = x
        y1 = y
    elif X['type_widge'].value==1:
        x1 = x
        y1 = -(y-y[0])
    elif X['type_widge'].value==2:
        x1 = -x
        y1 = -(y-y[0])/2
    
    #truncate data as required
    Bmin = X['Bmin_widge'].value
    Bmax = X['Bmax_widge'].value
    idx = (x1>=Bmin) & (x1<=Bmax) & (x1>0)
    y0 = y1[idx]*X['Cfactor']
    x0 = np.log10(x1[idx])
    
    Mn = y0[-1]
    yn = y0/Mn
    
    #get NUTS parameters
    nSample = X['nSample_widge'].value
    nTune = X['nTune_widge'].value
    
    print('Unmixing model is processing. Please wait .....')
    trace, _ = br.unmix(nB,x0,yn,nsample=nSample,tune=nTune)
    print('Unmixing model is complete')


    X['nB'] = nB
    X['x0'] = x0
    X['y0'] = y0
    X['yn'] = yn
    X['Mn'] = Mn
    X['trace'] = trace
    
    return X

def plot_model(X):
    
    nB = X['nB']
    x = X['x0']
    y = X['yn']
    Mn = X['Mn']
    trace = X['trace']
    ounits = X['ounits_widge'].value
    
    _ = br.plot_mixture(x,y,trace,nB,ounits,Mn=Mn)
   
    plt.tight_layout()
     
    outputfile = X['name_widge'].value+'_MODEL.pdf'
    plt.savefig(outputfile, dpi=300, bbox_inches="tight") 
    
    return X