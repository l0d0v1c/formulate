import pandas as pd
import re
import numpy as np

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
        






    

        




                        
        
        