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

Addition .add takes the name of the component, a raw formula and a dictionary of physical properties. The mixture is stored in .mixture
```
c.mixture
{'Water': {'H': 111.017, 'O': 55.508, 'Hf': -15865.97, 'rho': 1.0, 'N': 0},
 'Nitrogen': {'N': 71.394, 'Hf': 0.0, 'rho': 0.01, 'O': 0, 'H': 0},
 'Oxygen': {'O': 62.502, 'Hf': 0.0, 'rho': 0.01, 'N': 0, 'H': 0}}
 ```

 The atoms amount are converted in ATG/kg as well as physical properties (True|False) as per kg.

