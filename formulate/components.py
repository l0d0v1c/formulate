import pandas as pd
import re
import numpy as np
import pickle,base64,sys,io,gzip
import pysmiles

class components():
    def __init__(self,physical={'Hf':True},rond=3):
        """
            physical:   list of physical parameters additive (True) or not (False) 
                        or non additive and mass propotional (None),
                        default enthalpy of formation Hf (kJ/mol)
            rond : number of decimals in formula for readability
            Chemical element.
            

        """
        self.mixture=dict()
        self.atoms={'H': {'n': 1,
              'name': 'Hydrogen',
              'mass': 1.007941,
              'period': 1,
              'series': 1,
              'oxistates': '1*, -1'},
             'He': {'n': 2,
              'name': 'Helium',
              'mass': 4.002602,
              'period': 1,
              'series': 2,
              'oxistates': '*'},
             'Li': {'n': 3,
              'name': 'Lithium',
              'mass': 6.94,
              'period': 2,
              'series': 3,
              'oxistates': '1*'},
             'Be': {'n': 4,
              'name': 'Beryllium',
              'mass': 9.0121831,
              'period': 2,
              'series': 4,
              'oxistates': '2*'},
             'B': {'n': 5,
              'name': 'Boron',
              'mass': 10.811,
              'period': 2,
              'series': 5,
              'oxistates': '3*'},
             'C': {'n': 6,
              'name': 'Carbon',
              'mass': 12.01074,
              'period': 2,
              'series': 1,
              'oxistates': '4*, 2, -4*'},
             'N': {'n': 7,
              'name': 'Nitrogen',
              'mass': 14.006703,
              'period': 2,
              'series': 1,
              'oxistates': '5, 4, 3, 2, -3*'},
             'O': {'n': 8,
              'name': 'Oxygen',
              'mass': 15.999405,
              'period': 2,
              'series': 1,
              'oxistates': '-2*, -1'},
             'F': {'n': 9,
              'name': 'Fluorine',
              'mass': 18.998403163,
              'period': 2,
              'series': 6,
              'oxistates': '-1*'},
             'Ne': {'n': 10,
              'name': 'Neon',
              'mass': 20.1797,
              'period': 2,
              'series': 2,
              'oxistates': '*'},
             'Na': {'n': 11,
              'name': 'Sodium',
              'mass': 22.98976928,
              'period': 3,
              'series': 3,
              'oxistates': '1*'},
             'Mg': {'n': 12,
              'name': 'Magnesium',
              'mass': 24.3051,
              'period': 3,
              'series': 4,
              'oxistates': '2*'},
             'Al': {'n': 13,
              'name': 'Aluminium',
              'mass': 26.9815385,
              'period': 3,
              'series': 7,
              'oxistates': '3*'},
             'Si': {'n': 14,
              'name': 'Silicon',
              'mass': 28.0855,
              'period': 3,
              'series': 5,
              'oxistates': '4*, -4'},
             'P': {'n': 15,
              'name': 'Phosphorus',
              'mass': 30.973761998,
              'period': 3,
              'series': 1,
              'oxistates': '5*, 3, -3'},
             'S': {'n': 16,
              'name': 'Sulfur',
              'mass': 32.0648,
              'period': 3,
              'series': 1,
              'oxistates': '6*, 4, 2, -2'},
             'Cl': {'n': 17,
              'name': 'Chlorine',
              'mass': 35.4529,
              'period': 3,
              'series': 6,
              'oxistates': '7, 5, 3, 1, -1*'},
             'Ar': {'n': 18,
              'name': 'Argon',
              'mass': 39.948,
              'period': 3,
              'series': 2,
              'oxistates': '*'},
             'K': {'n': 19,
              'name': 'Potassium',
              'mass': 39.0983,
              'period': 4,
              'series': 3,
              'oxistates': '1*'},
             'Ca': {'n': 20,
              'name': 'Calcium',
              'mass': 40.078,
              'period': 4,
              'series': 4,
              'oxistates': '2*'},
             'Sc': {'n': 21,
              'name': 'Scandium',
              'mass': 44.955908,
              'period': 4,
              'series': 8,
              'oxistates': '3*'},
             'Ti': {'n': 22,
              'name': 'Titanium',
              'mass': 47.867,
              'period': 4,
              'series': 8,
              'oxistates': '4*, 3'},
             'V': {'n': 23,
              'name': 'Vanadium',
              'mass': 50.9415,
              'period': 4,
              'series': 8,
              'oxistates': '5*, 4, 3, 2, 0'},
             'Cr': {'n': 24,
              'name': 'Chromium',
              'mass': 51.9961,
              'period': 4,
              'series': 8,
              'oxistates': '6, 3*, 2, 0'},
             'Mn': {'n': 25,
              'name': 'Manganese',
              'mass': 54.938044,
              'period': 4,
              'series': 8,
              'oxistates': '7, 6, 4, 3, 2*, 0, -1'},
             'Fe': {'n': 26,
              'name': 'Iron',
              'mass': 55.845,
              'period': 4,
              'series': 8,
              'oxistates': '6, 3*, 2, 0, -2'},
             'Co': {'n': 27,
              'name': 'Cobalt',
              'mass': 58.933194,
              'period': 4,
              'series': 8,
              'oxistates': '3, 2*, 0, -1'},
             'Ni': {'n': 28,
              'name': 'Nickel',
              'mass': 58.6934,
              'period': 4,
              'series': 8,
              'oxistates': '3, 2*, 0'},
             'Cu': {'n': 29,
              'name': 'Copper',
              'mass': 63.546,
              'period': 4,
              'series': 8,
              'oxistates': '2*, 1'},
             'Zn': {'n': 30,
              'name': 'Zinc',
              'mass': 65.38,
              'period': 4,
              'series': 8,
              'oxistates': '2*'},
             'Ga': {'n': 31,
              'name': 'Gallium',
              'mass': 69.723,
              'period': 4,
              'series': 7,
              'oxistates': '3*'},
             'Ge': {'n': 32,
              'name': 'Germanium',
              'mass': 72.63,
              'period': 4,
              'series': 5,
              'oxistates': '4*'},
             'As': {'n': 33,
              'name': 'Arsenic',
              'mass': 74.921595,
              'period': 4,
              'series': 5,
              'oxistates': '5, 3*, -3'},
             'Se': {'n': 34,
              'name': 'Selenium',
              'mass': 78.971,
              'period': 4,
              'series': 1,
              'oxistates': '6, 4*, -2'},
             'Br': {'n': 35,
              'name': 'Bromine',
              'mass': 79.9035,
              'period': 4,
              'series': 6,
              'oxistates': '7, 5, 3, 1, -1*'},
             'Kr': {'n': 36,
              'name': 'Krypton',
              'mass': 83.798,
              'period': 4,
              'series': 2,
              'oxistates': '2*'},
             'Rb': {'n': 37,
              'name': 'Rubidium',
              'mass': 85.4678,
              'period': 5,
              'series': 3,
              'oxistates': '1*'},
             'Sr': {'n': 38,
              'name': 'Strontium',
              'mass': 87.62,
              'period': 5,
              'series': 4,
              'oxistates': '2*'},
             'Y': {'n': 39,
              'name': 'Yttrium',
              'mass': 88.90584,
              'period': 5,
              'series': 8,
              'oxistates': '3*'},
             'Zr': {'n': 40,
              'name': 'Zirconium',
              'mass': 91.224,
              'period': 5,
              'series': 8,
              'oxistates': '4*'},
             'Nb': {'n': 41,
              'name': 'Niobium',
              'mass': 92.90637,
              'period': 5,
              'series': 8,
              'oxistates': '5*, 3'},
             'Mo': {'n': 42,
              'name': 'Molybdenum',
              'mass': 95.95,
              'period': 5,
              'series': 8,
              'oxistates': '6*, 5, 4, 3, 2, 0'},
             'Tc': {'n': 43,
              'name': 'Technetium',
              'mass': 97.9072,
              'period': 5,
              'series': 8,
              'oxistates': '7*'},
             'Ru': {'n': 44,
              'name': 'Ruthenium',
              'mass': 101.07,
              'period': 5,
              'series': 8,
              'oxistates': '8, 6, 4*, 3*, 2, 0, -2'},
             'Rh': {'n': 45,
              'name': 'Rhodium',
              'mass': 102.9055,
              'period': 5,
              'series': 8,
              'oxistates': '5, 4, 3*, 1*, 2, 0'},
             'Pd': {'n': 46,
              'name': 'Palladium',
              'mass': 106.42,
              'period': 5,
              'series': 8,
              'oxistates': '4, 2*, 0'},
             'Ag': {'n': 47,
              'name': 'Silver',
              'mass': 107.8682,
              'period': 5,
              'series': 8,
              'oxistates': '2, 1*'},
             'Cd': {'n': 48,
              'name': 'Cadmium',
              'mass': 112.414,
              'period': 5,
              'series': 8,
              'oxistates': '2*'},
             'In': {'n': 49,
              'name': 'Indium',
              'mass': 114.818,
              'period': 5,
              'series': 7,
              'oxistates': '3*'},
             'Sn': {'n': 50,
              'name': 'Tin',
              'mass': 118.71,
              'period': 5,
              'series': 7,
              'oxistates': '4*, 2*'},
             'Sb': {'n': 51,
              'name': 'Antimony',
              'mass': 121.76,
              'period': 5,
              'series': 5,
              'oxistates': '5, 3*, -3'},
             'Te': {'n': 52,
              'name': 'Tellurium',
              'mass': 127.6,
              'period': 5,
              'series': 5,
              'oxistates': '6, 4*, -2'},
             'I': {'n': 53,
              'name': 'Iodine',
              'mass': 126.90447,
              'period': 5,
              'series': 6,
              'oxistates': '7, 5, 1, -1*'},
             'Xe': {'n': 54,
              'name': 'Xenon',
              'mass': 131.293,
              'period': 5,
              'series': 2,
              'oxistates': '2, 4, 6'},
             'Cs': {'n': 55,
              'name': 'Caesium',
              'mass': 132.90545196,
              'period': 6,
              'series': 3,
              'oxistates': '1*'},
             'Ba': {'n': 56,
              'name': 'Barium',
              'mass': 137.327,
              'period': 6,
              'series': 4,
              'oxistates': '2*'},
             'La': {'n': 57,
              'name': 'Lanthanum',
              'mass': 138.90547,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Ce': {'n': 58,
              'name': 'Cerium',
              'mass': 140.116,
              'period': 6,
              'series': 9,
              'oxistates': '4, 3*'},
             'Pr': {'n': 59,
              'name': 'Praseodymium',
              'mass': 140.90766,
              'period': 6,
              'series': 9,
              'oxistates': '4, 3*'},
             'Nd': {'n': 60,
              'name': 'Neodymium',
              'mass': 144.242,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Pm': {'n': 61,
              'name': 'Promethium',
              'mass': 144.9128,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Sm': {'n': 62,
              'name': 'Samarium',
              'mass': 150.36,
              'period': 6,
              'series': 9,
              'oxistates': '3*, 2'},
             'Eu': {'n': 63,
              'name': 'Europium',
              'mass': 151.964,
              'period': 6,
              'series': 9,
              'oxistates': '3*, 2'},
             'Gd': {'n': 64,
              'name': 'Gadolinium',
              'mass': 157.25,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Tb': {'n': 65,
              'name': 'Terbium',
              'mass': 158.92535,
              'period': 6,
              'series': 9,
              'oxistates': '4, 3*'},
             'Dy': {'n': 66,
              'name': 'Dysprosium',
              'mass': 162.5,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Ho': {'n': 67,
              'name': 'Holmium',
              'mass': 164.93033,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Er': {'n': 68,
              'name': 'Erbium',
              'mass': 167.259,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Tm': {'n': 69,
              'name': 'Thulium',
              'mass': 168.93422,
              'period': 6,
              'series': 9,
              'oxistates': '3*, 2'},
             'Yb': {'n': 70,
              'name': 'Ytterbium',
              'mass': 173.054,
              'period': 6,
              'series': 9,
              'oxistates': '3*, 2'},
             'Lu': {'n': 71,
              'name': 'Lutetium',
              'mass': 174.9668,
              'period': 6,
              'series': 9,
              'oxistates': '3*'},
             'Hf': {'n': 72,
              'name': 'Hafnium',
              'mass': 178.49,
              'period': 6,
              'series': 8,
              'oxistates': '4*'},
             'Ta': {'n': 73,
              'name': 'Tantalum',
              'mass': 180.94788,
              'period': 6,
              'series': 8,
              'oxistates': '5*'},
             'W': {'n': 74,
              'name': 'Tungsten',
              'mass': 183.84,
              'period': 6,
              'series': 8,
              'oxistates': '6*, 5, 4, 3, 2, 0'},
             'Re': {'n': 75,
              'name': 'Rhenium',
              'mass': 186.207,
              'period': 6,
              'series': 8,
              'oxistates': '7, 6, 4, 2, -1'},
             'Os': {'n': 76,
              'name': 'Osmium',
              'mass': 190.23,
              'period': 6,
              'series': 8,
              'oxistates': '8, 6, 4*, 3, 2, 0, -2'},
             'Ir': {'n': 77,
              'name': 'Iridium',
              'mass': 192.217,
              'period': 6,
              'series': 8,
              'oxistates': '6, 4*, 3, 2, 1*, 0, -1'},
             'Pt': {'n': 78,
              'name': 'Platinum',
              'mass': 195.084,
              'period': 6,
              'series': 8,
              'oxistates': '4*, 2*, 0'},
             'Au': {'n': 79,
              'name': 'Gold',
              'mass': 196.966569,
              'period': 6,
              'series': 8,
              'oxistates': '3*, 1'},
             'Hg': {'n': 80,
              'name': 'Mercury',
              'mass': 200.592,
              'period': 6,
              'series': 8,
              'oxistates': '2*, 1'},
             'Tl': {'n': 81,
              'name': 'Thallium',
              'mass': 204.3834,
              'period': 6,
              'series': 7,
              'oxistates': '3, 1*'},
             'Pb': {'n': 82,
              'name': 'Lead',
              'mass': 207.2,
              'period': 6,
              'series': 7,
              'oxistates': '4, 2*'},
             'Bi': {'n': 83,
              'name': 'Bismuth',
              'mass': 208.9804,
              'period': 6,
              'series': 7,
              'oxistates': '5, 3*'},
             'Po': {'n': 84,
              'name': 'Polonium',
              'mass': 208.9824,
              'period': 6,
              'series': 5,
              'oxistates': '6, 4*, 2'},
             'At': {'n': 85,
              'name': 'Astatine',
              'mass': 209.9871,
              'period': 6,
              'series': 6,
              'oxistates': '7, 5, 3, 1, -1*'},
             'Rn': {'n': 86,
              'name': 'Radon',
              'mass': 222.0176,
              'period': 6,
              'series': 2,
              'oxistates': '2*'},
             'Fr': {'n': 87,
              'name': 'Francium',
              'mass': 223.0197,
              'period': 7,
              'series': 3,
              'oxistates': '1*'},
             'Ra': {'n': 88,
              'name': 'Radium',
              'mass': 226.0254,
              'period': 7,
              'series': 4,
              'oxistates': '2*'},
             'Ac': {'n': 89,
              'name': 'Actinium',
              'mass': 227.0278,
              'period': 7,
              'series': 10,
              'oxistates': '3*'},
             'Th': {'n': 90,
              'name': 'Thorium',
              'mass': 232.0377,
              'period': 7,
              'series': 10,
              'oxistates': '4*'},
             'Pa': {'n': 91,
              'name': 'Protactinium',
              'mass': 231.03588,
              'period': 7,
              'series': 10,
              'oxistates': '5*, 4'},
             'U': {'n': 92,
              'name': 'Uranium',
              'mass': 238.02891,
              'period': 7,
              'series': 10,
              'oxistates': '6*, 5, 4, 3'},
             'Np': {'n': 93,
              'name': 'Neptunium',
              'mass': 237.0482,
              'period': 7,
              'series': 10,
              'oxistates': '6, 5*, 4, 3'},
             'Pu': {'n': 94,
              'name': 'Plutonium',
              'mass': 244.0642,
              'period': 7,
              'series': 10,
              'oxistates': '6, 5, 4*, 3'},
             'Am': {'n': 95,
              'name': 'Americium',
              'mass': 243.0614,
              'period': 7,
              'series': 10,
              'oxistates': '6, 5, 4, 3*'},
             'Cm': {'n': 96,
              'name': 'Curium',
              'mass': 247.0704,
              'period': 7,
              'series': 10,
              'oxistates': '4, 3*'},
             'Bk': {'n': 97,
              'name': 'Berkelium',
              'mass': 247.0703,
              'period': 7,
              'series': 10,
              'oxistates': '4, 3*'},
             'Cf': {'n': 98,
              'name': 'Californium',
              'mass': 251.0796,
              'period': 7,
              'series': 10,
              'oxistates': '4, 3*'},
             'Es': {'n': 99,
              'name': 'Einsteinium',
              'mass': 252.083,
              'period': 7,
              'series': 10,
              'oxistates': '3*'},
             'Fm': {'n': 100,
              'name': 'Fermium',
              'mass': 257.0951,
              'period': 7,
              'series': 10,
              'oxistates': '3*'},
             'Md': {'n': 101,
              'name': 'Mendelevium',
              'mass': 258.0984,
              'period': 7,
              'series': 10,
              'oxistates': '3*'},
             'No': {'n': 102,
              'name': 'Nobelium',
              'mass': 259.101,
              'period': 7,
              'series': 10,
              'oxistates': '3, 2*'},
             'Lr': {'n': 103,
              'name': 'Lawrencium',
              'mass': 262.1096,
              'period': 7,
              'series': 10,
              'oxistates': '3*'},
             'Rf': {'n': 104,
              'name': 'Rutherfordium',
              'mass': 267.1218,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
             'Db': {'n': 105,
              'name': 'Dubnium',
              'mass': 268.1257,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
             'Sg': {'n': 106,
              'name': 'Seaborgium',
              'mass': 271.1339,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
             'Bh': {'n': 107,
              'name': 'Bohrium',
              'mass': 272.1383,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
             'Hs': {'n': 108,
              'name': 'Hassium',
              'mass': 270.1343,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
             'Mt': {'n': 109,
              'name': 'Meitnerium',
              'mass': 276.1516,
              'period': 7,
              'series': 8,
              'oxistates': '*'},
              
              }

        self.physical=physical
        self.rond=rond
        self.allatoms=set()
        
    def add(self,name,rawformula,param):
        """
        add a component having rawformula and param, convert it is ATG/kg
        
        """
        if rawformula!='':

          r=re.findall('([A-Z][a-z]?)([0-9\.]+)?',rawformula)
          molecule={g[0]:(1 if g[1]=='' else float(g[1])) for g in r}
          MassMolar=lambda molec:sum([self.atoms[at]['mass']*n for at,n in molec.items()])
          MM=MassMolar(molecule)
          for at,n in molecule.items():
              self.allatoms.add(at)
              molecule[at]=round(n*1000/MM,self.rond) 
        else:
          molecule={}

        for p,addit in self.physical.items():
            if p not in param.keys():
                raise Exception("The physical value for %s is required" %(p,))

            molecule[p]=round(param[p]*1000/MM if addit!=None else param[p],self.rond) 
        self.mixture.update({name:molecule})
        # Update mixture missing atoms
        for name,m in self.mixture.items():
            for at in self.allatoms:
                if at not in m.keys():
                    m[at]=0
            self.mixture[name]=m

    def setrates(self,rates):
        """
        Define the amount of each component eg. rates={"water":1.0}
        """
        total=0
        self.rates=dict()
        for comp,rate in rates.items():
            if comp not in self.mixture.keys():
                raise Exception("The component %s is not in the mixture" %(comp,))
            total+=rate
        for comp,rate in rates.items():
            self.rates[comp]=round(rate/total,self.rond)

    def mixing(self):
        """
        return formulation
        self.formulationtab : pandas
        self.formulation: dict
        """
        tab,col=[],['Component','Rate']
        for name,_ in self.mixture.items():
            l=[name,self.rates[name]]
            for at in self.allatoms:
                if at not in col:
                    col.append(at)
                l.append(self.mixture[name][at])

            for param,addit in self.physical.items():
                if param not in col:
                    col.append(param)
                l.append(self.mixture[name][param])
            
            tab.append(l)
            col2=['Component','Rate']
            for i in range(2,len(col)):
                cl=col[i]

                col2.append(cl)
        
        tab=pd.DataFrame(tab,columns=col2)
        tabm=['Formulation',1.0]

        self.formulation=dict()
        for col in tab.columns[2:]:
            if col in self.allatoms or self.physical[col]:
                s=0
                for x,rate in zip(tab[col],tab['Rate']):
                    s+=x*rate
                tabm.append(s)
                self.formulation[col]=s
            else:
                tabm.append("Non additive")
                self.formulation[col]="NA"

        self.formulationtab=pd.DataFrame([tabm],columns=tab.columns)

        self.formulationlist=pd.concat([tab,self.formulationtab])
        self.formulationlist.reset_index(inplace=True,drop=True)
        

        return self.formulationlist

    def oxygenbalance(self):
        """
        Oxygen balance of the formulation for CHON mixtures in %
        """
        if len(self.allatoms-set(['C','H','O','N']))>0:
            raise Exception("Only valid for CHON molecules")
        return round(-1600/1000*(2*self.formulation['C']+self.formulation['H']/2-self.formulation['O']),self.rond)

    def eutectic(self,underrelax=0.01):
        """
          param:underrelax : underrelaxation for Newton-Raphson method
          The formulation must have Hfus in J/mol and Tfus in K for each component
          The row formula is not required

          Brunet, L., J. Caillard, et P. André. « Thermodynamic calculation of n-component 
          eutectic mixtures ». International Journal of Modern Physics C 15, nᵒ 5 (2004): 675‑87. 
          https://doi.org/10.1142/S0129183104006121.



        """
        nc=len(self.mixture)
        T=self.formulationlist["Tfus"].values[:-1].mean()
     
        x=[1/nc for _ in range(nc-1)]
        x.append(T)
        GO=True
        
        while GO:
            mat=[]
            for compo in range(nc-1):
                l=[0 for _ in range(nc)]
                l[compo]=1/x[compo]

                l[-1]=-self.formulationlist["Hfus"].values[compo]/8.314/T/T
                mat.append(l)
            l=[-1/(1-sum(x[:-1])) for _ in range(nc-1)]
            l.append(-self.formulationlist["Hfus"].values[-2]/8.314/T/T)
            mat.append(l)
            mat=np.linalg.inv(np.array(mat))
            l=[[np.log(x[compo])+self.formulationlist["Hfus"].values[compo]/8.314/T
                -self.formulationlist["Hfus"].values[compo]/8.314/self.formulationlist["Tfus"].values[compo]] for compo in range(nc-1)]
            l.append([np.log(1-sum(x[:-1]))+self.formulationlist["Hfus"].values[-2]/8.314/T
                -self.formulationlist["Hfus"].values[-2]/8.314/self.formulationlist["Tfus"].values[-2]])
            l=np.array(l)
            dx=np.matmul(mat,l)
            x=[x[i]-underrelax*dx[i,0] for i in range(nc)]
            T=x[-1]
            s=0
            for xx in dx:
                s+=xx*xx
            GO= False if s<0.1 else True
        x=x[:-1]
        x.append(1-sum(x))
        x.append(1)
        self.formulationlist["Rate"]=pd.DataFrame([[round(xx,self.rond)] for xx in x])
        TT=self.formulationlist["Tfus"].values
        TT[-1]=round(T,0)
        self.formulationlist["Tfus"]=TT


        return round(T,0)
    
    def equilibrium(self):
      """
          Calculation of temperature of equilibrium by an experimental feedforward neural network
          only for C, H, O, N molecules
          Return temperature in K



      """
      if len(self.allatoms-set(['C','H','O','N']))>0:
        raise Exception("Only valid for CHON molecules")
      neuro='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'
      modele=pickle.load(io.BytesIO(gzip.decompress(base64.b64decode(neuro.encode()))))
      iw,ow,biases=modele[4],modele[5],modele[6]
      boundaries={'C':[1,77],'H':[0,100],'O':[0,61],'N':[0,70.57],'Hf':[-499993,499993],'Tf':[162.49,4687.95]}
      mixt=self.formulation
      mixt["Hf"]/=4.18
      eq= (lambda x,b: ((b[0]+b[1])/2-b[0])*x+(b[0]+b[1])/2)((lambda X:np.dot((lambda X:(lambda x: 1. / (1 + np.exp(-x)))(np.dot(X, iw)+biases))([X[j] for j in ['C',"H","O","N","Hf"]]),ow))((lambda mixture,b: {at:(lambda x,b:(x-(b[0]+b[1])/2)/((b[0]+b[1])/2-b[0]))(v,b[at]) for at,v in mixt.items()})(mixt,boundaries)),boundaries['Tf'])
      return round(eq,0)


    def enthalpyf(self,smiles,hamiltonian="PM3"):
      """
          Enthalpy of formation at 298K gas phase by hamiltoniam PM3 and AM1 (Mopac) in kJ/mol for
          C,H,O,N molecules
          param:smiles: Canonical smiles molecule encoding
          param:hamiltonian: PM3 or AM1

      """
      if len(self.allatoms-set(['C','H','O','N']))>0:
        raise Exception("Only valid for CHON molecules")

      npm3='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'
      nam1='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'
      modele=pickle.load(io.BytesIO(gzip.decompress(base64.b64decode(npm3.encode() if hamiltonian=="PM3" else nam1.encode()))))
      iw,ow,biases=modele[4],modele[5],modele[6]
      mol_with_H = pysmiles.read_smiles(smiles, explicit_hydrogen=True)
      v=[0 for _ in range(len(modele[2])-1)]
      for c in mol_with_H.nodes(data='element'):
          v[list(modele[2]).index(c[1])]+=1
      for symb in modele[2][4:-1]:
          v[list(modele[2]).index(symb)]=smiles.count(symb)
      
      my_dict = {k: v for k, v in modele[0].items() }
      boundaries={i:[j[1],j[2]] for i,j in my_dict.items()} 
      #print(boundaries)
      #print(v)
      recode=lambda x,b:(x-(b[0]+b[1])/2)/((b[0]+b[1])/2-b[0])
      for i in range(0,len(v)):
          v[i]=recode(v[i],boundaries[modele[2][i]])
      #print(v)
      predict=lambda X:np.dot((lambda X:(lambda x: 1. / (1 + np.exp(-x)))(np.dot(X, iw)+biases))(X),ow)
      hr=predict(v)
      #print(hr)
      decode=lambda x,b: ((b[0]+b[1])/2-b[0])*x+(b[0]+b[1])/2
      
      
      return decode(hr,boundaries["Hf%s" % hamiltonian])*4.18/1000


        
        
        






    

        




                        
        
        