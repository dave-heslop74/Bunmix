import ipysheet
import ipywidgets as ipyw
import matplotlib.pyplot as plt
from ipysheet import column, row
import numpy as np
import matplotlib.ticker as mtick

def enter_data(X):
    nval = 500
    x = np.zeros(nval)
    x[:] = np.nan
    y = np.zeros(nval)
    y[:] = np.nan

    sheet= ipysheet.sheet(rows=nval,columns=2,column_headers=(('Field [mT]','Magnetization')))
    sheet.layout.height = '300px'
    sheet.layout.width = '600px'

    col1 = column(0,x)
    col2 = column(1,y)
    
    display(sheet)
    
    X['sheet']=sheet
    
    return X

def load_data(X):
    df = ipysheet.to_dataframe(X['sheet'])
    B = np.array(df['Field [mT]'])
    B = B[~np.isnan(B)]

    M = np.array(df['Magnetization'])
    M = M[~np.isnan(M)]
    
    X['B'] = B
    X['M'] = M
    
    return X

def get_data(X):
    style = {'description_width': 'initial'} #general style settings

    
    X = load_data(X)
    x = np.log10(X['B'])
    y = X['M']
    
    X['name_widge'] = ipyw.Text(style={'description_width': 'initial'},
    value='my sample',
    description='Sample name:',
    disabled=False)

    X['mass_widge'] = ipyw.FloatText(style={'description_width': 'initial'},
    value='-1.0',
    description='Sample mass [g]:',
    disabled=False)

    #SAMPLE VOLUME widget
    X['volume_widge'] = ipyw.FloatText(style={'description_width': 'initial'},
    value='-1.0',
    description='Sample volume [cm$^3$]:',
    disabled=False)
    
    # INPUT UNIT widget
    X['iunits_widge'] = ipyw.Dropdown(style={'description_width': 'initial'},
    options=[('emu',0),('emu/cm'+u'\u00B3',1),('emu/g',2),('Am'+u'\u00B2', 3),('A/m', 4),('Am'+u'\u00B2'+'/kg', 5)],
    value=0,
    description='Input Magnetization Units:')
    
    # OUTPUT UNIT widget
    X['ounits_widge'] = ipyw.Dropdown(style={'description_width': 'initial'},
    options=[('emu',0),('emu/cm'+u'\u00B3',1),('emu/g',2),('Am'+u'\u00B2', 3),('A/m', 4),('Am'+u'\u00B2'+'/kg', 5)],
    value=0,
    description='Output Magnetization Units:')

    # SHOW PLOT widget
    X['show_widge'] = ipyw.Checkbox(value=False, description = 'Plot data')
    
    # SAVE PLOT widget
    X['save_widge'] = ipyw.Checkbox(value=False, description = 'Save plot') 

    # MIN & MAX FIELD widget
    X['Bmin_widge'] = ipyw.FloatText(value=np.min(10**x),description='Minimum B [mT]',step=1,style=style)    
    X['Bmax_widge'] = ipyw.FloatText(value=np.max(10**x),description='Maximum B [mT]',step=1,style=style)
    
    p = ipyw.interactive(plot_data,
                         x=ipyw.fixed(x),
                         y=ipyw.fixed(y),
                         show = X['show_widge'],
                         name = X['name_widge'],
                         mass = X['mass_widge'],
                         volume = X['volume_widge'],
                         ounits=X['ounits_widge'],
                        iunits=X['iunits_widge'],
                        Bmin=X['Bmin_widge'],
                         Bmax=X['Bmax_widge'],
                         save=X['save_widge'])

    tab_nest = ipyw.Tab()
    tab_nest.set_title(0, 'Input data')
    tab_nest.children = [ipyw.VBox(children = p.children)]
    display(tab_nest)
    
    #display(p)
    
    return X
    
def plot_data(x,y,show,name,mass,volume,iunits,ounits,Bmin,Bmax,save):
    
    idx = (10**x>=Bmin) & (10**x<=Bmax)
    y0 = y[idx]
    x0 = x[idx]
    
    error = 0
    Cfactor = 1
    
    ##### UNIT CONVERSION FROM EMU #####
    # Convert emu to emu/cm3
    if (iunits==0) & (ounits==1):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = 1/volume
    
    
    # Convert emu to emu/g
    if (iunits==0) & (ounits==2):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = 1/mass
    
    # Convert emu to Am2
    if (iunits==0) & (ounits==3):
            error = 0
            Cfactor = 0.001    
    
    # Convert emu to A/m
    if (iunits==0) & (ounits==4):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = 0.001/(volume/1E6)
            
    # Convert emu to Am2/kg
    if (iunits==0) & (ounits==5):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = 1/mass

    ##### UNIT CONVERSION FROM EMU/CM3 #####
    # Convert emu/cm3 to emu
    if (iunits==1) & (ounits==0):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = volume
    
    
    # Convert emu/cm3 to emu/g
    if (iunits==1) & (ounits==2):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = volume/mass
    
    # Convert emu/cm3 to Am2
    if (iunits==1) & (ounits==3):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = volume/0.001    
    
    # Convert emu/cm3 to A/m
    if (iunits==1) & (ounits==4):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = volume*0.001
            
    # Convert emu/cm3 to Am2/kg
    if (iunits==1) & (ounits==5):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = volume/mass

    ##### UNIT CONVERSION FROM EMU/G ###############################
    # Convert emu/g to emu
    if (iunits==2) & (ounits==0):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = mass
    
    # Convert emu/g to emu/cm3
    if (iunits==2) & (ounits==1):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = mass/volume
    
    # Convert emu/g to Am2
    if (iunits==2) & (ounits==3):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = mass*0.001    
    
    # Convert emu/g to A/m
    if (iunits==2) & (ounits==4):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
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
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = 1000/volume
    
    # Convert Am2 to emu/g
    if (iunits==3) & (ounits==2):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = 1000/mass    
    
    # Convert Am2 to A/m
    if (iunits==3) & (ounits==4):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = 1/(volume/1E6)
            
    # Convert Am2 to Am2/kg
    if (iunits==3) & (ounits==5):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = 1/(mass/1000)

    ##### UNIT CONVERSION FROM A/m ###############################
    # Convert A/m to emu
    if (iunits==4) & (ounits==0):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = 1000*volume
    
    # Convert A/m to emu/cm3
    if (iunits==4) & (ounits==1):
            error = 0
            Cfactor = 1000
    
    # Convert A/m to emu/g
    if (iunits==4) & (ounits==2):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = 1000*volume/mass    
    
    # Convert A/m to Am2
    if (iunits==4) & (ounits==3):
        if volume<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample volume</h2>"))
        else:
            error = 0
            Cfactor = volume/1E6
            
    # Convert A/m to Am2/kg
    if (iunits==4) & (ounits==5):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = (volume/1E6)/(mass/1000)

    ##### UNIT CONVERSION FROM Am2/kg ###############################
    # Convert Am2/kg to emu
    if (iunits==5) & (ounits==0):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = mass
    
    # Convert Am2/kg to emu/cm3
    if (iunits==5) & (ounits==1):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = mass/volume
    
    # Convert Am2/kg to emu/g
    if (iunits==5) & (ounits==2):
        if mass<=0:
            error = 0
            Cfactor = 1    
    
    # Convert Am2/kg to Am2
    if (iunits==5) & (ounits==3):
        if mass<=0:
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass</h2>"))
        else:
            error = 0
            Cfactor = mass/1E3
            
    # Convert Am2/kg to A/m
    if (iunits==5) & (ounits==4):
        if (mass<=0) | (volume<=0):
            error = 1
            display(ipyw.HTML(value="<h2>Error: Conversion requires a sample mass and volume</h2>"))
        else:
            error = 0
            Cfactor = (mass/1000)/(volume/1E6)
            
    ##### FINISHED UNIT CONVERSION ######
    #X['Cfactor'] = Cfactor
    if error == 0:
        fig0, (ax1, ax2) = plt.subplots(1,2,figsize=([13, 4.8]))
        ax1.plot(10**x0,y0*Cfactor,'.k-')
        ax1.set_xscale('log')
        ax1.set_xlim([np.min(10**x[idx]),np.max(10**x[idx])])
        ax1.set_xlabel('B [mT]',fontsize=14)
    
        if ounits==0:
            ax1.set_ylabel('Magnetic moment [emu]',fontsize=14)
        elif ounits==1:
            ax1.set_ylabel('Magnetization [emu / cm'+u'\u00B3'+']',fontsize=14)
        elif ounits==2:
            ax1.set_ylabel('Mass magnetization [emu / g]',fontsize=14)
        elif ounits==3:
            ax1.set_ylabel('Magnetic moment [Am'+u'\u00B2'+']',fontsize=14)  
        elif ounits==4:
            ax1.set_ylabel('Magnetization [A / m]',fontsize=14)
        else:
            ax1.set_ylabel('Mass magnetization [Am'+u'\u00B2'+' / kg]',fontsize=14)
        
        ax1.tick_params(labelsize=14)
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        ylim = ax1.get_ylim()
        if np.min(y0*Cfactor)>=0.0:
            ax1.set_ylim([0,ylim[1]])
        
        
        Dy = np.diff(y0*Cfactor)/np.diff(x0)
        Dx = x0[:-1]+np.diff(x0)/2
        ax2.plot(10**Dx,Dy,'.k-')
        ax2.set_xscale('log')
        ax2.set_xlim([np.min(10**x[idx]),np.max(10**x[idx])])
        ax2.set_xlabel('B [mT]',fontsize=14)
        
        if ounits==0:
            ax2.set_ylabel('Derivative [emu / $\log_{10}$(B[mT])]',fontsize=14)
        elif ounits==1:
            ax2.set_ylabel('Derivative [emu / cm'+u'\u00B3'+' / $\log_{10}$(B[mT])]',fontsize=14)
        elif ounits==2:
            ax2.set_ylabel('Derivative [emu / g / $\log_{10}$(B[mT])]',fontsize=14)
        elif ounits==3:
            ax2.set_ylabel('Derivative [Am'+u'\u00B2'+' / $\log_{10}$(B[mT])]',fontsize=14)  
        elif ounits==4:
            ax2.set_ylabel('Derivative [A / m / $\log_{10}$(B[mT])]',fontsize=14)
        else:
            ax2.set_ylabel('Derivative [Am'+u'\u00B2'+' / kg / $\log_{10}$(B[mT])]',fontsize=14)
    
        ax2.tick_params(labelsize=14)
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    
        plt.tight_layout()
    
        if save==True:
            outputfile = name+'_DATA.pdf'
            plt.savefig(outputfile, dpi=300, bbox_inches="tight")       
        
    return 1