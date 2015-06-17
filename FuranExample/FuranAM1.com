#P AM1(Input,print) CIS(Singlets,AllTransitionDensities,NStates=15) pop(full) force

Thiophene AM1 CIS energy

0   1
6        0.955475   -0.715963    0.000000
6       -0.345627   -1.092818    0.000000
8       -1.157305    0.000000    0.000000
6       -0.345627    1.092818    0.000000
6        0.955475    0.715963    0.000000
1        1.810438    1.369382    0.000000
1       -0.840307    2.046884    0.000000
1       -0.840307   -2.046884    0.000000
1        1.810438   -1.369382    0.000000

Method=8 CoreType=1 PeptideFC=0.0114505045 RIJScale=0.5291772086
****
H
PQN=1 NValence=1 F0ss=0.5725046400 ZetaOverlap=1.5161372400 U=-0.3898159000
Beta=-0.1468379300 CoreKO=0.9367259700 KON=0,0,0,1.0220707700 KON=1,0,1,0.8573527300
KON=0,1,1,1.1167017200 KON=2,1,1,0.9565689600 EISol=-0.4188109900 EHeat=0.0830298200
Alpha=1.5322399100
GCore=0.0072864600,1.4810771900,2.8898973500
GCore=0.0003116200,1.1395427000,2.8617012500
GCore=-0.0011016100,0.4660914200,2.3989195800
****
C
PQN=2,2 NValence=4 F0ss=0.4182154200 F0sp=0.4544492400 F0pp=0.3452489400 F2pp=0.1901027700
G1sp=0.2635705200
ZetaOverlap=2.1213905800,1.6393033300
U=-1.8354384900,-1.3854253300
Beta=-0.5700795100,-0.3627690700 DDN=0,1,0.7960081100 DDN=1,1,0.7320042300
CoreKO=1.1394499500 KON=0,0,0,0.7141943600 KON=1,0,1,0.7517453400 KON=0,1,1,1.1767038100
KON=2,1,1,0.8159232100 EISol=-4.4398989700 EHeat=0.2723305500 Alpha=1.3934555100
DipHyp=1.5919246299
GCore=0.0008540200,1.3325270400,3.4177634700
GCore=0.0032292100,1.2262351500,3.8210561300
GCore=-0.0014695500,1.2404556300,5.2374863600
GCore=-0.0000813100,1.2708306300,5.9841122600
****
O
PQN=2,2 NValence=6 F0ss=0.6120795600 F0sp=0.5321302400 F0pp=0.4942333100 F2pp=0.2673750200
G1sp=0.4626171900
ZetaOverlap=3.2677171200,2.7061243100
U=-3.6308722500,-2.9511601900
Beta=-1.1878632300,-1.1537714000 DDN=0,1,0.4985610200 DDN=1,1,0.4526235800
CoreKO=0.9614567800 KON=0,0,0,0.8312087100 KON=1,0,1,0.5027290600 KON=0,1,1,0.7782572500
KON=2,1,1,0.6076144000 EISol=-11.6164445900 EHeat=0.0949133100 Alpha=2.5956664100
DipHyp=0.9970646747
GCore=0.0209894200,0.8811997200,1.8384538900
GCore=0.0065002600,1.8521608100,3.1242482700
****

