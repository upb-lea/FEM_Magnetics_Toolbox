Include "PreParameter.pro";


//Rx = Rc/Sin[Pi/3]/Fill; // original wrong formula
Rx = Rc*Sqrt[2*Pi/3/Sqrt[3]/Fill]; // corrected by Till

Rx = Rx; // work in meters not millimeters anymore
Ry = Rx;
Rc  = Rc; // work in meters not millimeters anymore

AreaCond = Pi*Rc^2;
AreaCell = AreaCond/Fill;


NbrCond = 3*(NbrLayers+1)^2 - 3*(NbrLayers+1) + 1 ; // only hexagonal case

DefineConstant[
  MD = {1/1.5, Name "Cell/99Mesh density factor"} // original definition
];

MD = 1; // to make it fast !attention! high means fast/coarse


p = Rc/10*MD  ;
pc = Rc/10*MD ;


ISOL  =  2000;
COND  = 10000;
BOUND = 20000;