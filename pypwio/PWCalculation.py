import pprint

from .PWOutput import PWOutput
from .PWXML import PWXML
from .PWInput import PWInput

class PWCalculation():
    '''
    Data and information necessary to run a pw.x calculation through nanoHUB
    :param pwinput: PWInput object from pymatgen
    :param pwoutput: PWOutput object from htdft
    :param metadata: information for nanoHUB `submit` and i/o parsing
        required values ['input_path', 'output_path', 'pseudodir_path']
        values with defaults {'walltime': '01:00:00', 'ncpus': 1}
        optional values ['run_name', 'espresso_name']
    '''
    def __init__(self, pwinput, metadata, pwoutput=None, pwxml=None):
        self.pwinput = pwinput
        self.pwoutput = pwoutput
        self.pwxml = pwxml
        self.metadata = metadata

        return

    def __repr__(self):
        return 'METADATA:\n{}\n\nINPUT:\n"""\n{}\n"""\n\nOUTPUT:\n{}\n\nXML:\n{}\n'.format(pprint.pformat(self.metadata),
                                                                                str(self.pwinput),
                                                                                self.pwoutput,
                                                                                self.pwxml)

    def write_input(self, path=None):
        if path:
            self.metadata['input_path'] =  path
        elif 'input_path' in self.metadata.keys():
            if self.metadata['input_path']:
                path = self.metadata['input_path']
            else:
                raise ValueError('Must provide a valid filepath in function argument or metadata[\'input_path\']')

        try:
            self.pwinput.write_file(self.metadata['input_path'])
            return self.metadata['input_path']
        except Exception as e:
            print(e)
            return

    def parse_output(self, path=None):
        if path:
            self.metadata['output_path'] = path
        elif 'output_path' in self.metadata:
            if self.metadata['output_path']:
                path = self.metadata['output_path']
            else:
                raise ValueError('Must provide a valid filepath in function argument or metadata[\'output_path\']')
        with open(path, 'r') as f:
            text = f.read()
        self.pwoutput = PWOutput(stdout_string=text, stdout_path=path)
        if self.pwoutput:
            self.metadata['DFTman_status'].append('parsed')
        return self.pwoutput

    def parse_xml(self, path=None):
        if path:
            self.metadata['xml_path'] = path
        elif 'xml_path' in self.metadata:
            if self.metadata['xml_path']:
                path = self.metadata['xml_path']
            else:
                raise ValueError('Must provide a valid filepath in function argument or metadata[\'xml_path\']')
        with open(path, 'r') as f:
            text = f.read()
        self.pwxml = PWXML(xml_string=text, xml_path=path)
        return self.pwxml

    def as_dict(self):
        if self.pwoutput:
            output = self.pwoutput.as_dict()
        else:
            output = None
        if self.pwxml:
            xml = self.pwxml.as_dict()
        else:
            xml = None
        pwdict = {
            'input': self.pwinput.as_dict(),
            'output': output,
            'xml': xml,
            'metadata': self.metadata
        }
        return pwdict

    @classmethod
    def from_dict(cls, pwdict):
        pwinput = PWInput.from_dict(pwdict['input'])
        if pwdict['output']:
            pwoutput = PWOutput.from_dict(pwdict['output'])
        else:
            pwoutput = pwdict['output']
        if pwdict['xml']:
            pwxml = PWXML.from_dict(pwdict['xml'])
        else:
            pwxml = pwdict['xml']

        metadata = pwdict['metadata']
        pwcalculation = cls(pwinput, metadata, pwoutput, pwxml)
        return pwcalculation
