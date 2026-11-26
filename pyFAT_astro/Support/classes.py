


class Galaxy:
    def __init__(self, name = 'Empty'):
        self.name = 'Empty'
        self.center = []
        self.Tirific_Template = None
        self.config = None
        self.fits_files = None
        self.fable_dict = None
        self.hdr = None
        self.wcs = None
        self.filename = None
        self.directory = None
        self.channel_width = None



    def get_model_parameters(self,variables = None):
        if variables is None:
            variables = []
        model_parameters = {}
        for var in variables:
            model_parameters[var] = getattr(self, var, None)
        return model_parameters

