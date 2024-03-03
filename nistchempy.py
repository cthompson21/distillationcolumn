'''
Python API for NIST Chemistry WebBook
'''

#%% Imports

import sys
if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources

import re, os, requests, zipfile
from urllib.parse import urlparse, parse_qs
from bs4 import BeautifulSoup, Comment
import pandas as pd
import pdb


#%% Package info

__version__ = '0.2.3'




#%% All NIST Data
ASCII_BETA = '\u03B2'
ASCII_ALPHA  = '\u03B1'


def get_all_data():
    '''
    Returns pandas dataframe containing info on all NIST Chem WebBook compounds
    '''
    dt0 = {'mol_weight': 'float64'}
    dt1 = {k: 'string' for k in ('ID', 'name', 'formula', 'inchi', 'inchi_key', 'cas_rn')}
    dt2 = {k: 'bool' for k in ('mol2D', 'mol3D', 'cIR', 'cTZ', 'cMS', 'cUV', 'cGC',
                               'cTG', 'cTC', 'cTP', 'cSO', 'cTR', 'cIE', 'cIC', 'cES', 'cDI')}
    dtypes = {**dt0, **dt1, **dt2}
    pkg = importlib_resources.files('nistchempy')
    data_file = pkg / 'nist_data.zip'
    with importlib_resources.as_file(data_file) as path:
        zf = zipfile.ZipFile(path) 
        df = pd.read_csv(zf.open('nist_data.csv'), dtype = dtypes)
        zf.close()
    
    return df


#%% Support functions

def _is_compound(soup):
    '''
    Checks if html is a single compound page and returns NIST ID if yes
    '''
    header = soup.findAll('h1', {'id': 'Top'})
    if not header:
        return None
    # get info
    header = header[0]
    info = header.findNext('ul')
    if not info:
        return None
    # extract NIST ID
    for comment in soup.findAll(text = lambda text: isinstance(text, Comment)):
        comment = str(comment).replace('\r\n', '').replace('\n', '')
        if not '/cgi/cbook.cgi' in comment:
            continue
        return re.search(r'/cgi/cbook.cgi\?Form=(.*?)&', comment).group(1)
    
    return None


#%% Compound-related classes

class Spectrum():
    '''
    Class for IR, MS, and UV-Vis extracted from NIST Chemistry WebBook
    '''
    
    _pretty_names = {'IR': 'IR spectrum',
                     'TZ': 'THz IR spectrum',
                     'MS': 'Mass spectrum',
                     'UV': 'UV-Vis spectrum'}
    
    def __init__(self, compound, spec_type, spec_idx, jdx):
        self.compound = compound
        self.spec_type = spec_type
        self.spec_idx = spec_idx
        self.jdx_text = jdx
    
    def save(self, name = None, path_dir = None):
        '''
        Saves spectrum in JDX format
        '''
        path = name if name else f'{self.compound.ID}_{self.spec_type}_{self.spec_idx}.jdx'
        if path_dir:
            path = os.path.join(path_dir, path)
        with open(path, 'w') as outf:
            outf.write(self.jdx_text)
    
    def __str__(self):
        return f'Spectrum({self.compound.ID}, {self._pretty_names[self.spec_type]} #{self.spec_idx})'
    
    def __repr__(self):
        return f'Spectrum({self.compound.ID}, {self._pretty_names[self.spec_type]} #{self.spec_idx})'


class Compound():
    '''
    Object for NIST Chemistry WebBook compound
    '''
    
    # NIST URLs
    _NIST_URL = 'https://webbook.nist.gov'
    _COMP_ID = '/cgi/cbook.cgi'
    
    # mappings for spectra
    _MASKS = {'1': 'cTG', '2': 'cTC', '4': 'cTP', '8': 'cTR', '10': 'cSO',
              '20': 'cIE', '40': 'cIC', '80': 'cIR', '100': 'cTZ', '200': 'cMS',
              '400': 'cUV', '800': 'cES', '1000': 'cDI', '2000': 'cGC'}
    _SPECS = {'IR': 'IR', 'TZ': 'THz', 'MS': 'Mass', 'UV': 'UVVis'}
    
    def _load_compound_info(self):
        '''
        Loads main compound info
        '''
        # pdb.set_trace()
        if self.ID!=None:
            request_input = self._NIST_URL + self._COMP_ID, {'ID': self.ID, 'Units': 'SI'}
        elif self.name != None:
            request_input = self._NIST_URL + self._COMP_ID, {'Name': self.name,'Units': 'SI'}
        r = requests.get(*request_input)
        
        
        if not r.ok:
            raise ConnectionError(f'Bad NIST response, status code: {r.status_code}')
        # check if it is compound page
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
        # pdb.set_trace()
        header = soup.findAll('h1', {'id': 'Top'})
        
        if not header:
            raise ValueError(f'Bad compound ID: {self.ID}')
        header = header[0]
        # get info
        info = header.findNext('ul')
        if not info:
            raise ValueError(f'Bad compound ID: {self.ID}')
        # name
        self.name = header.text
        # synonyms
        hits = info.findChildren(text = re.compile('Other names'))
        if hits:
            text = hits[0].findParent('li').text.replace('Other names:', '')
            synonyms = [_.strip(';').strip() for _ in text.split('\n')]
            self.synonyms = [_ for _ in synonyms if _]
        # formula
        hits = info.findChildren(text = re.compile('Formula'))
        if hits:
            text = hits[0].findParent('li').text.replace('Formula:', '')
            self.formula = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', text.strip())
        # mol weight
        hits = info.findChildren(text = re.compile('Molecular weight'))
        if hits:
            text = hits[0].findParent('li').text.replace('Molecular weight:', '')
            self.mol_weight = float(text)
        # InChI and InChI key
        hits = info.findChildren(attrs = {'class': 'inchi-text'})
        if hits:
            for hit in hits:
                if 'InChI=' in hit.text:
                    self.inchi = hit.text
                elif re.search(r'', hit.text):
                    self.inchi_key = hit.text
        # CAS RN
        hits = info.findChildren(text = re.compile('CAS Registry Number:'))
        if hits:
            text = hits[0].findParent('li').text.replace('CAS Registry Number:', '')
            self.cas_rn = text.strip()
        # 2D structure
        hits = info.findChildren(attrs = {'href': re.compile('Str2File')})
        if hits:
            self.data_refs['mol2D'] = self._NIST_URL + hits[0].attrs['href']
        # 3D structure
        hits = info.findChildren(attrs = {'href': re.compile('Str3File')})
        if hits:
            self.data_refs['mol3D'] = self._NIST_URL + hits[0].attrs['href']
        # other data and spectroscopy
        hits = info.findChildren(attrs = {'href': re.compile('/cgi/cbook.cgi.*Mask=\d')})
        for hit in hits:
            mask = re.search('Mask=(\d+)', hit.attrs['href']).group(1)
            key = self._MASKS.get(mask, hit.text)
            if key in self.data_refs:
                self.data_refs[key] += [self._NIST_URL + hit.attrs['href']]
            else:
                self.data_refs[key] = [self._NIST_URL + hit.attrs['href']]
    
    def get_2D(self):
        '''
        Loads 2D structure in MOL2 format
        '''
        if 'mol2D' not in self.data_refs:
            return
        r = requests.get(self.data_refs['mol2D'])
        if r.ok:
            self.mol2D = r.text
    
    def get_3D(self):
        '''
        Loads 3D structure in MOL2 format
        '''
        if 'mol3D' not in self.data_refs:
            return
        r = requests.get(self.data_refs['mol3D'])
        if r.ok:
            self.mol3D = r.text
    
    def get_spectra(self, spec_type):
        '''
        Loads available mass spectra in JCAMP-DX format
        '''
        
        
        if spec_type not in self._SPECS:
            raise ValueError(f'Bad spec_type value: {spec_type}')
        if 'c'+spec_type not in self.data_refs:
            return
        r = requests.get(self.data_refs['c'+spec_type][0])
        if not r.ok:
            return
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
        # get available spectrum indexes
        idxs = soup.findAll(attrs = {'href': re.compile('Index=')})
        idxs = [re.search(r'Index=(\d+)', _.attrs['href']).group(1) for _ in idxs]
        idxs = sorted(list(set(idxs)))
        # load jdxs
        for idx in idxs:
            spec = requests.get(self._NIST_URL + self._COMP_ID,
                                {'JCAMP': self.ID, 'Index': idx,
                                 'Type': self._SPECS[spec_type]})
            if spec.ok:
                spec = Spectrum(self, spec_type, idx, spec.text)
                getattr(self, spec_type).append(spec)
                
    def get_thermo_condensed(self):
        '''
        Loads thermochemistry data
        '''
        spec_type='IR'
        
       
        if spec_type not in self._SPECS:
            raise ValueError(f'Bad spec_type value: {spec_type}')
        # if 'c'+spec_type not in self.data_refs:
        #     return
        
        
        if self.ID!=None:
            request_input = self._NIST_URL + self._COMP_ID, {'ID': self.ID, 'Units': 'SI', 'Mask':'2#Thermo-Condensed'}
        elif self.name != None:
            request_input = self._NIST_URL + self._COMP_ID, {'Name': self.name,'Units': 'SI', 'Mask':'2#Thermo-Condensed'}
        r = requests.get(*request_input)
        
        # r = requests.get("https://webbook.nist.gov/cgi/cbook.cgi?ID=C110543&Units=SI&Mask=2#Thermo-Condensed")
        if not r.ok:
            return
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
        
        # get available spectrum indexes

        idxs = soup.findAll(attrs={'href':re.compile('Thermo-Phase')})
        for idx in idxs[0:1]:

            r = requests.get(self._NIST_URL + idx.attrs['href'])
        if not r.ok:
            return 
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
            
        idxs = soup.findAll('table', )
        titles=[]
        
        #### make dictionary with the table data and table title
        for idx in idxs:
            
            if 'aria-label' in idx.attrs.keys():
                to_append = idx.attrs['aria-label']
                while  to_append in titles:
                    to_append+='_1'
                    
                titles.append(to_append)
            else:
                titles.append('')
        tables=pd.read_html(r.text)
        
        tables = dict(zip(titles, tables))
        
        
        ##### here, get the heat of vaporization
        tit=None
        ev_table=None
        single_temp=False
        import operator
        temp ='\t'.join(tables.keys())
        if 'Enthalpy of vaporization_1' in temp:
            tit='Enthalpy of vaporization_1'
            ev_table=tables[tit]
        elif 'Enthalpy of vaporization' in temp:   
            tit='Enthalpy of vaporization'
            ev_table=tables[tit]
        else:
            ev_table==None
        
        

        if 'Tc (K)' in np.array(ev_table) or 'Tc (K)' in ev_table.columns  :
                single_temp=False
        elif "vapH (kJ/mol)" in '\t'.join(ev_table.columns):
                single_temp=True
                
                self.delta_H_vap_300K = ev_table.to_numpy()[0,0]
                self.delta_vap_H_constants = pd.DataFrame(columns = ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc'],
                                                          data = [[self.delta_H_vap_300K,-1,-1,-1]],
                                                          index=[0])
        else:
            print("Cant decipher the tables")
        
            
        if single_temp==False:    
                
                
            if "Temperature (K)" in np.array(ev_table):
                ev_table.index = ev_table[0]
                
                ev_table = ev_table.drop(columns=[0])
                ev_table = ev_table.rename(mapper = {1:'value'},axis=1)
            elif "Temperature (K)" in ev_table.columns:
                
                ev_table = ev_table.transpose()
                ev_table = ev_table.rename(mapper = {0:'value'},axis=1)
                
            
            x = re.findall('\d+[\s-]*\d+', ev_table['value']["Temperature (K)"])
            min_temp=int(x[0])
            max_temp = int(x[1])
            
            A = float(ev_table['value']['A (kJ/mol)'])
            beta = float(ev_table['value'][ASCII_BETA])
            if ASCII_ALPHA in ev_table.index:
                alpha = float(ev_table['value'][ASCII_BETA])
            else:
                alpha=beta
    
            Tc=float(ev_table['value']['Tc (K)'])
            T=np.arange(min_temp,max_temp+1)
            Tr = T/Tc
            self.delta_H_vap_300K = A*np.exp(-alpha*300/Tc)*(1-300/Tc)**beta
            self.delta_vap_H_constants = pd.DataFrame(columns = ['Enthalpy Vaporization (kJ/mol)','delta_H_vap_A','delta_H_vap_beta','delta_H_vap_Tc'],
                                                      data = [[delta_H_vap_300K,A,beta,Tc]],
                                                      index=[0])
            
                    
            self.delta_vap_H = pd.DataFrame(data=A*np.exp(-alpha*Tr)*(1-Tr)**beta, index = T)
           
            
        ###now, get the vapor pressures.  
        tit=None
        self.p_vap = pd.DataFrame(columns = ['p_vap','T','Tmin','Tmax', 'ref'])
        if 'Antoine Equation Parameters' in tables.keys():
            
            if 'Temperature (K)' in tables['Antoine Equation Parameters'].columns:
                
                #log10(P) = A − (B / (T + C))
                
                
                for row in range(len(tables['Antoine Equation Parameters'])):
                    try:
                        x = re.findall('[\d.][\s-]*[\d.]+', tables['Antoine Equation Parameters'].iloc[row]["Temperature (K)"])
                        if not x:
                            continue
                        if isinstance(x,str):
                            print('x is a string')
                        if type(x) == list:
                            if len(x)!=2:
                                print('x is a list with length <2', x)
                            
                    except:
                        raise
                            
                    # print(tables['Antoine Equation Parameters'].iloc[row]["Temperature (K)"],x)
                    min_temp=float(x[0])
                    max_temp = float(x[1])
                    T=np.arange(min_temp, max_temp, 0.01)
                    A = tables['Antoine Equation Parameters'].iloc[row]['A']
                    
                    B=tables['Antoine Equation Parameters'].iloc[row]['B']
                    C=tables['Antoine Equation Parameters'].iloc[row]['C']
                    p=10**(A-(B/(T+C)))
                    ref = np.array([tables['Antoine Equation Parameters'].iloc[row]["Reference"]]*len(p))
                   
                    min_temp=np.array([min_temp,]*len(p))
                    max_temp=np.array([max_temp,]*len(p))
                    
                    self.antoine_parameters = [A,B,C]
                    
                    data = np.transpose([p,T,min_temp,max_temp, ref]) 
                    self.p_vap = pd.concat((self.p_vap, pd.DataFrame(index = T, data = data, columns = self.p_vap.columns)), axis = 0)
                    
            # pdb.set_trace()
            #### round the temperature values to two digits and eliminate the duplicates
            self.p_vap['T'] = self.p_vap['T'].apply(float).round(1)
            self.p_vap = self.p_vap.drop_duplicates(subset=['T'])
            
            
            
            ## reindex by temperature now that duplicate temperatures were eliminated
            self.p_vap = self.p_vap.loc[self.p_vap['p_vap'].isnull()==False]
            
            self.p_vap.index = self.p_vap['T']
            
            #pdb.set_trace()
            
            new_T = np.round(np.arange(np.min(self.p_vap.index), np.max(self.p_vap.index), 0.1),1)
            self.p_vap = self.p_vap[['p_vap', 'ref']]
            p_vap2 = pd.DataFrame(index = new_T, )
            p_vap2 = pd.concat((p_vap2, self.p_vap), axis = 1)# = self.p_vap[self.p_vap2.index]
            self.p_vap = p_vap2.copy()
            self.p_vap['p_vap']=self.p_vap['p_vap'].apply(float)
            self.p_vap=self.p_vap.sort_index()
            
            self.p_vap['ref'].loc[self.p_vap['p_vap'].isnull()]='interpolated'
            
            
            self.p_vap['p_vap'] =self.p_vap['p_vap'].interpolate(method = 'index')
            
            
            
            
            
            #### now get heat capacities
    def get_heat_capacities(self):
        
         
        if self.ID!=None:
            request_input = self._NIST_URL + self._COMP_ID, {'ID': self.ID, 'Units': 'SI', 'Mask':'2#Thermo-Condensed'}
        elif self.name != None:
            request_input = self._NIST_URL + self._COMP_ID, {'Name': self.name,'Units': 'SI', 'Mask':'2#Thermo-Condensed'}
        r = requests.get(*request_input)
        if not r.ok:
            return
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
            
        idxs = soup.findAll('table', )
        titles=[]
        
        #### make dictionary with the table data and table title
        for idx in idxs:
            
            if 'aria-label' in idx.attrs.keys():
                to_append = idx.attrs['aria-label']
                while  to_append in titles:
                    to_append+='_1'
                    
                titles.append(to_append)
            else:
                titles.append('')
        tables=pd.read_html(r.text)
        
        tables = dict(zip(titles, tables))
        
        tit=None
      
        
        if 'Constant pressure heat capacity of liquid' in tables.keys():
            
            if 'Temperature (K)' in tables['Constant pressure heat capacity of liquid'].columns:
                print(tables['Constant pressure heat capacity of liquid'])
                self.Cp_l = tables['Constant pressure heat capacity of liquid']['Cp,liquid (J/mol*K)']
                print("Found Cp_l", self.Cp_l)
         
                

     
        if self.ID!=None:
            request_input = self._NIST_URL + self._COMP_ID, {'ID': self.ID, 'Units': 'SI', 'Mask':'1#Thermo-Gas'}
        elif self.name != None:
            request_input = self._NIST_URL + self._COMP_ID, {'Name': self.name,'Units': 'SI', 'Mask':'1#Thermo-Gas'}
        r = requests.get(*request_input)
        if not r.ok:
            return
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
            
        idxs = soup.findAll('table', )
        titles=[]
        
        #### make dictionary with the table data and table title
        for idx in idxs:
            
            if 'aria-label' in idx.attrs.keys():
                to_append = idx.attrs['aria-label']
                while  to_append in titles:
                    to_append+='_1'
                    
                titles.append(to_append)
            else:
                titles.append('')
        tables=pd.read_html(r.text)
        
        tables = dict(zip(titles, tables))
        
        tit=None
      
        
        if 'Constant pressure heat capacity of gas' in tables.keys():
            
            if 'Temperature (K)' in tables['Constant pressure heat capacity of gas'].columns:
                print(tables['Constant pressure heat capacity of gas'].columns)
                self.Cp_g = tables['Constant pressure heat capacity of gas']['Cp,gas (J/mol*K)']
                print("Found Cp_g", self.Cp_g)
        else:
            print("No tables appropriate for CpG found")
         
                
                   
            
        return 
        
        def get_ir_spectra(self):
            '''
            Loads available IR spectra in JCAMP-DX format
            '''
            
            return self.get_spectra('IR')
    
    def get_tz_spectra(self):
        '''
        Loads available IR spectra in JCAMP-DX format
        '''
        
        return self.get_spectra('TZ')
    
    def get_ms_spectra(self):
        '''
        Loads available mass spectra in JCAMP-DX format
        '''
        
        return self.get_spectra('MS')
    
    def get_uv_spectra(self):
        '''
        Loads available UV-Vis spectra in JCAMP-DX format
        '''
        
        return self.get_spectra('UV')
    
    def get_all_spectra(self):
        '''
        Loads available spectroscopic data
        '''
        self.get_ir_spectra()
        self.get_tz_spectra()
        self.get_ms_spectra()
        self.get_uv_spectra()
    
    def get_all_data(self):
        '''
        Loads available structural and spectroscopic data
        '''
        self.get_2D()
        self.get_3D()
        self.get_all_spectra()
    
    def save_spectra(self, spec_type, path_dir = './'):
        '''
        Saves all spectra of given type to the specified folder
        '''
        if not os.path.isdir(path_dir):
            raise ValueError(f'"{path_dir}" must be directory')
        for spec in getattr(self, spec_type):
            spec.save(f'{self.ID}_{spec_type}_{spec.spec_idx}.jdx', path_dir)
    
    def save_ir_spectra(self, path_dir = './'):
        '''
        Saves IR spectra to the specified folder
        '''
        self.save_spectra('IR', path_dir)
    
    def save_tz_spectra(self, path_dir = './'):
        '''
        Saves IR spectra to the specified folder
        '''
        self.save_spectra('TZ', path_dir)
    
    def save_ms_spectra(self, path_dir = './'):
        '''
        Saves mass spectra to the specified folder
        '''
        self.save_spectra('MS', path_dir)
    
    def save_uv_spectra(self, path_dir = './'):
        '''
        Saves all UV-Vis spectra to the specified folder
        '''
        self.save_spectra('UV', path_dir)
    
    def save_all_spectra(self, path_dir = './'):
        '''
        Saves all UV-Vis spectra to the specified folder
        '''
        self.save_ir_spectra(path_dir)
        self.save_tz_spectra(path_dir)
        self.save_ms_spectra(path_dir)
        self.save_uv_spectra(path_dir)
    
    def __init__(self, search_value, search_option='ID'):
        if search_option=='ID':
            self.ID = search_value
        
            for prop, val in [('name', None), ('synonyms', []), ('formula', None), ('mol_weight', None),
                          ('inchi', None), ('inchi_key', None), ('cas_rn', None),
                          ('IR', []), ('TZ', []), ('MS', []), ('UV', []),
                          ('mol2D', None), ('mol3D', None),
                          ('data_refs', {})]:
                setattr(self, prop, val)
        elif search_option=='NAME':
            self.ID=None
            for prop, val in [('name', search_value), ('synonyms', []), ('formula', None), ('mol_weight', None),
                          ('inchi', None), ('inchi_key', None), ('cas_rn', None),
                          ('IR', []), ('TZ', []), ('MS', []), ('UV', []),
                          ('mol2D', None), ('mol3D', None),
                          ('data_refs', {})]:
                setattr(self, prop, val)
            
        self._load_compound_info()
    
    def __str__(self):
        return f'Compound({self.ID})'
    
    def __repr__(self):
        return f'Compound({self.ID})'


#%% Search-related classes

def print_search_parameters():
    '''
    Prints available search parameters
    '''
    info = {'Units': 'Units for thermodynamic data, "SI" or "CAL" for calorie-based',
            'MatchIso': 'Exactly match the specified isotopes (formula search only)',
            'AllowOther': 'Allow elements not specified in formula (formula search only)',
            'AllowExtra': 'Allow more atoms of elements in formula than specified (formula search only)',
            'NoIon': 'Exclude ions from the search (formula search only)',
            'cTG': 'Contains gas-phase thermodynamic data',
            'cTC': 'Contains condensed-phase thermodynamic data',
            'cTP': 'Contains phase-change thermodynamic data',
            'cTR': 'Contains reaction thermodynamic data',
            'cIE': 'Contains ion energetics thermodynamic data',
            'cIC': 'Contains ion cluster thermodynamic data',
            'cIR': 'Contains IR data',
            'cTZ': 'Contains THz IR data',
            'cMS': 'Contains MS data',
            'cUV': 'Contains UV/Vis data',
            'cGC': 'Contains gas chromatography data',
            'cES': 'Contains vibrational and electronic energy levels',
            'cDI': 'Contains constants of diatomic molecules',
            'cSO': 'Contains info on Henry\'s law'}
    max_len = max([len(_) for _ in info])
    spaces = [' '*(max_len - len(_) + 1) for _ in info]
    for (key, val), space in zip(info.items(), spaces):
        print(f'{key}{space}:   {val}')


class SearchParameters():
    '''
    Object containing parameters for compound search in NIST WebBook
    To get the full description of available options please use
    "print_search_parameters" function
    '''
    info = {'Units': 'SI',
            'MatchIso': False, 'AllowOther': False, 'AllowExtra': False, 'NoIon': False,
            'cTG': False, 'cTC': False, 'cTP': False, 'cTR': False, 'cIE': False, 'cIC': False,
            'cIR': False, 'cTZ': False, 'cMS': False, 'cUV': False, 'cGC': False,
            'cES': False, 'cDI': False, 'cSO': False}
    
    def get_request_parameters(self):
        '''
        Returns dictionary with GET parameters
        '''
        params = {'Units': self.Units}
        for key in self.info:
            if key == 'Units':
                continue
            val = getattr(self, key)
            if val:
                params[key] = 'on'
        
        return params
    
    def __init__(self, **kwargs):
        # set default
        for key, val in self.info.items():
            setattr(self, key, val)
        # check kwargs
        for key, val in kwargs.items():
            if key not in self.info:
                raise TypeError(f'"{key}" is an invalid keyword argument for SearchParameters')
            if key == 'Units' and val not in ('SI', 'CAL'):
                raise ValueError(f'Bad value for "Units" parameter: {val}')
            if key != 'Units' and type(val) is not bool:
                raise ValueError(f'Bad value for "{key}" parameter: {val}')
            setattr(self, key, val)
    
    def __str__(self):
        sep = ', ' # ',\n' + ' '*17
        text = [f'SearchParameters(Units={self.Units}'] + \
               [f'{key}={getattr(self, key)}' for key in self.info if key != 'Units' and getattr(self, key)]
        text[-1] = text[-1] + ')'
        
        return sep.join(text)
    
    def __repr__(self):
        sep = ', ' # ',\n' + ' '*17
        text = [f'SearchParameters(Units={self.Units}'] + \
               [f'{key}={getattr(self, key)}' for key in self.info if key != 'Units' and getattr(self, key)]
        text[-1] = text[-1] + ')'
        
        return sep.join(text)


class Search():
    '''
    Object for searching in NIST Chemistry WebBook
    '''
    
    # NIST URLs
    _NIST_URL = 'https://webbook.nist.gov'
    _COMP_ID = '/cgi/cbook.cgi'
    
    # parameters data
    search_types = {'formula': 'Formula', 'name': 'Name', 'inchi': 'InChI', 'cas': 'ID'}
    formula_only = ('MatchIso', 'AllowOther', 'AllowExtra', 'NoIon')
    
    def __init__(self, **kwargs):
        self.parameters = SearchParameters(**kwargs)
        self.IDs = []
        self.compounds = []
        self.lost = False
        self.success = True
    
    def find_compounds(self, identifier, search_type, **kwargs):
        '''
        Search for species data by chemical name
        search_type must be one of 'formula', 'name', 'inchi', 'cas'
        clear_found: clear all found compounds
        raise_lost: raise exception if limit of 400 compounds per search
                    was achieved
        '''
        if search_type not in self.search_types:
            raise ValueError(f'Bad search_type value: {search_type}')
        # prepare GET parameters
        params = {self.search_types[search_type]: identifier}
        params.update(self.parameters.get_request_parameters())
        addend = SearchParameters(**kwargs)
        params.update(addend.get_request_parameters())
        # load webpage
        r = requests.get(self._NIST_URL + self._COMP_ID, params)
        if not r.ok:
            self.success = False
            self.IDs = []
            self.compounds = []
            self.lost = False
            return
        soup = BeautifulSoup(re.sub('clss=', 'class=', r.text),
                             features = 'html.parser')
        # check if no compounds
        if search_type == 'inchi':
            errs = ['information from the inchi', 'no matching species found']
        else:
            errs = ['not found']
        err_flag = False
        for err in errs:
            if sum([err in _.text.lower() for _ in soup.findAll('h1')]):
                err_flag = True
                break
        if err_flag:
            self.success = True
            self.IDs = []
            self.compounds = []
            self.lost = False
            return
        # check if one compound
        flag = _is_compound(soup)
        if flag:
            self.success = True
            self.IDs = [flag]
            self.compounds = []
            self.lost = False
            return
        # extract IDs
        refs = soup.find('ol').findChildren('a', href = re.compile(self._COMP_ID))
        IDs = [parse_qs(urlparse(a.attrs['href']).query)['ID'][0] for a in refs]
        self.IDs = IDs
        self.compounds = []
        self.success = True
        self.lost = 'Due to the large number of matching species' in soup.text
    
    def load_found_compounds(self):
        '''
        Loads compounds
        '''
        self.compounds = [Compound(ID) for ID in self.IDs]
    
    def __str__(self):
        return f'Search(Success={self.success}, Lost={self.lost}, Found={len(self.IDs)})'
    
    def __repr__(self):
        return f'Search(Success={self.success}, Lost={self.lost}, Found={len(self.IDs)})'


import numpy as np
from matplotlib.pyplot import *
# X=Compound('C110543')  ##hexane
if __name__ == '__main__':
    Hvap = pd.DataFrame(index = np.arange(0,1000), columns = ['α-Methylstyrene', 'styrene'])
        # 'methane','ethane',
        #                                                       'propane','butane', 'pentane','hexane',
        #                                                       'heptane','octane',
        #                                                       'methanol', 'ethanol', 'n-propanol', 'n-butanol',
        #                                                       'toluene', 'acetic acid','sulfuric acid','α-Methylstyrene'])
    
    
    fig,ax = subplots(1,1)
    for alkane in Hvap.columns:
        X=Compound(alkane, search_option='NAME')  ##hexane
        print(X.name)
        X.get_heat_capacities()
        
        tables=X.get_thermo_condensed()
        print(X.p_vap['p_vap'])
        X.p_vap['p_vap']=X.p_vap['p_vap'].apply(float)
        X.p_vap['p_vap'].plot()
        suptitle('vapor pressure')
        
        
    
        
    legend( Hvap.columns)

        