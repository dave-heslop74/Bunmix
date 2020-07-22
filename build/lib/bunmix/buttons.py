import ipywidgets as ipyw
from traitlets import traitlets
import bunmix.input as bi
import bunmix.model as bm

class LoadedButton(ipyw.Button):
    """A button that can holds a value as a attribute."""

    def __init__(self, value=None, *args, **kwargs):
        super(LoadedButton, self).__init__(*args, **kwargs)
        # Create the value attribute.
        self.add_traits(value=traitlets.Any(value))
    
def data(ex):    
    bi.get_data(ex.value)
    
def params(ex):    
    bm.get_params(ex.value)
    
def model(ex):    
    bm.get_model(ex.value)
    
def plot(ex):    
    bm.plot_model(ex.value)