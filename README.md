# Formulate

Formulate is a library to build and manipulate formulations. It can be use for materials, cosmetics or any activities
involving mixing of components. This version computes oxygen balance and eutectic points.

## Basic object

The main object is 'components' and called by from components import components 
```
c=components(physical={"Hf":True,"rho":None}) 
c.add("Water","H2O",{'Hf':-285.83,"rho":1.0}) 
c.add("Nitrogen","N2",{'Hf':0,"rho":0.01}) 
c.add("Oxygen","O2",{'Hf':0,"rho":0.01})
```

The physical dictionary contains a list of physical properties that can be set True for the one that are additives (like enthaply, False when they depends on mass but are not additives (like heat of fusion), None when they are mass independant and non-additive (like temperature)

Addition .add takes the name of the component, a raw formula and a dictionary of physical properties. The raw formula may include all possible atoms in any order. The mixture is stored in .mixture
```
c.mixture
{'Water': {'H': 111.017, 'O': 55.508, 'Hf': -15865.97, 'rho': 1.0, 'N': 0},
 'Nitrogen': {'N': 71.394, 'Hf': 0.0, 'rho': 0.01, 'O': 0, 'H': 0},
 'Oxygen': {'O': 62.502, 'Hf': 0.0, 'rho': 0.01, 'N': 0, 'H': 0}}
 ```

 The atoms amount are converted in ATG/kg as well as physical properties (True|False) as per kg.
 The rates are defined by .setrates
 ```
 c.setrates({"Water":0.5,"Nitrogen":0.25,"Oxygen":0.25})
 c.mixing()
 ```
 .mixing outputs a table of the formulation

 ## Processed formula

 The following processed data are available:
 * c.formulationtab : a pandas dataframe of the final formulation
 * c.formulationlist : a pandas dataframe of components
 * c.formulation : a formulation dictionary

 ## Models

 ### Oxygen balance 
 c.oxygenbalance() gives the oxygen balance for C,H,N,O mixtures (https://en.wikipedia.org/wiki/Oxygen_balance)

 ### Eutectic
c.eutectic(underrelax=0.01) allow the computation of eutectic points. In this case, raw formula is not required but the physical properties must include Tfus (melting point in K) and Hfus (melting enthalpy in J/mol). For instance to find the eutectic point of LiF/NaF/KF mixture:
```
eutec=components(physical={"Hfus":None,"Tfus":None})
eutec.add("KF","",{"Hfus":28500,"Tfus":856+273})
eutec.add("LiF","",{"Hfus":10000,"Tfus":845+273})
eutec.add("NaF","",{"Hfus":32600,"Tfus":985+273})
eutec.setrates({"KF":1,"LiF":1,"NaF":1})
eutec.mixing()
T=eutec.eutectic()
eutec.formulationlist
```

underrelax defines the under relaxation of Newton-Raphson algorithm. Please refere to:
Brunet, L., J. Caillard, et P. André. « Thermodynamic calculation of n-component 
eutectic mixtures ». International Journal of Modern Physics C 15, nᵒ 5 (2004): 675‑87. 
https://doi.org/10.1142/S0129183104006121.



