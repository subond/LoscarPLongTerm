        LOSCAR MODEL


Loscar.m         script, main
LoscarDif.m      function, differential equations (DEQs)

 
THmfun.m         function, conveyor circulation
dafunPE.m        function, CO2 routine
kspfun.m         function, CaCO3 solubilities
PEmisfun.m       function, CO2 emissions (time interpolation)


./myode          solver DEQ, Matlab built-in functions
./dat            model input/output


./docs

Zeebe and Zachos, PA 2007
Zachos et al., Nature 2008
Zeebe et al., Science 2008
Zeebe et al., Science 2008 Supplement
Walker and Kasting, P^3, 1992
Toggweiler 1999

SchemeLoscar.jpg  schematic figure of Loscar model 
toggpo4.m         script, Toggweiler PO4 3-box model


%======= Check values

Zeebe et al. (2008) 
Supplement, Fig S1, green lines

Total emissions = 1000 Pg C
Release time    =  500 years

pCO2Max_pHAInit_pHAMin_DeltapH =
0.5050    8.1646    7.9561    0.2085


