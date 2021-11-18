//Include "cell_dat.pro";

//p  = Rc/3;
//pc = Rc/3;

Dx = 0.; Dy = 0.;

Geometry.AutoCoherence = 0;
Mesh.Algorithm = 1;


//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
Function HexCirc_
  dP=newp-1;
  dR=news-1;

  Point(dP+1)  = {Dx+Rx*Cos[1*Pi/6],  Dy+Ry*Sin[1*Pi/6], 0 , p};
  Point(dP+2)  = {Dx+Rx*Cos[3*Pi/6],  Dy+Ry*Sin[3*Pi/6], 0 , p};
  Point(dP+3)  = {Dx+Rx*Cos[5*Pi/6],  Dy+Ry*Sin[5*Pi/6], 0 , p};
  Point(dP+4)  = {Dx+Rx*Cos[7*Pi/6],  Dy+Ry*Sin[7*Pi/6], 0 , p};
  Point(dP+5)  = {Dx+Rx*Cos[9*Pi/6],  Dy+Ry*Sin[9*Pi/6], 0 , p};
  Point(dP+6)  = {Dx+Rx*Cos[11*Pi/6], Dy+Ry*Sin[11*Pi/6], 0 , p};
  Point(dP+7)  = {Dx, Dy, 0 , p};

  Line(dR+1) = {dP+1,dP+2};
  Line(dR+2) = {dP+2,dP+3};
  Line(dR+3) = {dP+3,dP+4};
  Line(dR+4) = {dP+4,dP+5};
  Line(dR+5) = {dP+5,dP+6};
  Line(dR+6) = {dP+6,dP+1};

  Line Loop(dR+1) = {dR+6,dR+1,dR+2,dR+3,dR+4,dR+5};

  Point(dP+8)  = {Dx, Dy, 0 , pc};
  Point(dP+9)  = {Dx+Rc, Dy, 0 , pc};
  Point(dP+10) = {Dx-Rc, Dy, 0 , pc};
  Point(dP+11) = {Dx, Dy+Rc, 0 , pc};
  Point(dP+12) = {Dx, Dy-Rc, 0 , pc};

  Circle(dR+7)  = {dP+9,dP+8,dP+11};
  Circle(dR+8)  = {dP+11,dP+8,dP+10};
  Circle(dR+9)  = {dP+10,dP+8,dP+12};
  Circle(dR+10) = {dP+12,dP+8,dP+9};

  Line Loop(dR+2) = {dR+7,dR+8,dR+9,dR+10};

  Plane Surface(dR+3) = {dR+1,dR+2};
  Plane Surface(dR+4) = {dR+2};

  surfCond[] += {dR+4};
  surfIsol[] += {dR+3};
Return

//=========================================

surfCond[] = {};
surfIsol[] = {};

Call HexCirc_ ;

If (NbrLayers > 0)
  dx1 = 2*Rx*Cos[Pi/6] ;
  dy1 = 2*Ry*Cos[Pi/6] ;
  iL = 0 ; kk = 0 ; l = 0 ;

  For iL In {1:NbrLayers}
    Dx = -iL * dx1 ;
    Dy = 0;
    For k In {1:6}
      For l In {1:iL}
        kk = 0 ;
        If (iL == NbrLayers)
          kk = k ;
        EndIf
        // COND += 1 ;
        // ISOL += 1 ;
        Dx += dx1 * Cos[(2-k)*Pi/3] ;
        Dy += dy1 * Sin[(2-k)*Pi/3] ;
        Call HexCirc_ ;
      EndFor
    EndFor
  EndFor
EndIf


Geometry.AutoCoherence = 1;
Coherence;

For k In {0:#surfCond[]-1}
  Physical Surface(COND+k) = {surfCond[k]};
  Physical Surface(ISOL+k) = {surfIsol[k]};
EndFor


allSurfaces[] = Surface '*' ;
lines_boundary[] = CombinedBoundary{Surface{allSurfaces[]}; };
Physical Line(BOUND) = lines_boundary[];
