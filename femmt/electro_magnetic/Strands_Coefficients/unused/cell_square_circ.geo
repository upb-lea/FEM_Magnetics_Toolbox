//Include "cell_square_circ_dat.pro";

cen0 = newp ; Point(cen0)  = { 0,  0, 0 , pc};

pw[]+=newp ; Point(pw[0])  = { Rc,  0, 0 , pc};
pw[]+=newp ; Point(pw[1])  = {  0, Rc, 0 , pc};
pw[]+=newp ; Point(pw[2])  = {-Rc,  0, 0 , pc};
pw[]+=newp ; Point(pw[3])  = {  0,-Rc, 0 , pc};

cw[]+=newl ; Circle(cw[0])  = {pw[0], cen0, pw[1]};
cw[]+=newl ; Circle(cw[1])  = {pw[1], cen0, pw[2]};
cw[]+=newl ; Circle(cw[2])  = {pw[2], cen0, pw[3]};
cw[]+=newl ; Circle(cw[3])  = {pw[3], cen0, pw[0]};

llcw[] += newll ; Line Loop(llcw[0]) = {cw[]};
sw[] += news ; Plane Surface(sw[0]) = {llcw[0]};

pb[]+=newp ; Point(pb[0])  = { Dex/2,  Dey/2, 0 , p};
pb[]+=newp ; Point(pb[1])  = {-Dex/2,  Dey/2, 0 , p};
pb[]+=newp ; Point(pb[2])  = {-Dex/2, -Dey/2, 0 , p};
pb[]+=newp ; Point(pb[3])  = { Dex/2, -Dey/2, 0 , p};

lb[]+=newl ; Line(lb[0])  = {pb[0], pb[1]};
lb[]+=newl ; Line(lb[1])  = {pb[1], pb[2]};
lb[]+=newl ; Line(lb[2])  = {pb[2], pb[3]};
lb[]+=newl ; Line(lb[3])  = {pb[3], pb[0]};

llb[] += newll ; Line Loop(llb[0]) = {lb[]};
is[] += news ; Plane Surface(is[0]) = {llb[0],llcw[0]};

If(NbrLayers==1)
  xaux[] = {0,      0, Dex, Dex,  Dex, -Dex, -Dex, -Dex};
  yaux[] = {Dey, -Dey, Dey,   0, -Dey,  Dey,    0, -Dey};

  For k In {0:7}
    surf[]=Translate {xaux[k], yaux[k], 0} {Duplicata{ Surface{sw[0]}; }}; sw[] += surf[0] ;
    surf[]=Translate {xaux[k], yaux[k], 0} {Duplicata{ Surface{is[0]}; }}; is[] += surf[0] ;
  EndFor
EndIf

For k In {0:#sw[]-1}
  Physical Surface(COND+k) = {sw[k]};
  Physical Surface(ISOL+k) = {is[k]};
EndFor

bndi[] =CombinedBoundary{Surface{is[],sw[]};};
Physical Line(BOUND) = {bndi[]};
