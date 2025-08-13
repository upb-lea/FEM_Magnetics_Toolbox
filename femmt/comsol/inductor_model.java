/*
 * inductor_model.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 25 2024, 14:47 by COMSOL 6.1.0.357. */
public class inductor_model {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("C:\\Users\\tpiepe\\Repositories\\FEM_Magnetics_Toolbox\\femmt\\comsol");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 2);
    model.component("comp1").geom("geom1").axisymmetric(true);

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").physics().create("mef", "ElectricInductionCurrents", "geom1");

    model.study().create("std1");
    model.study("std1").create("freq", "Frequency");
    model.study("std1").feature("freq").setSolveFor("/physics/mef", true);

    model.component("comp1").geom("geom1").create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("pol1").set("source", "table");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -0.001}, {0.002, -0.001}, {0.002, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.002, 0.002}, {0.002, 0.001}, {0, 0.001}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("pol1");
    model.component("comp1").geom("geom1").create("c1", "Circle");
    model.component("comp1").geom("geom1").feature("c1").set("r", 3.0E-4);
    model.component("comp1").geom("geom1").feature("c1").set("pos", new double[]{0.0025, 0});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -0.001}, {0.002, -0.001}, {0.002, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.002, 0.002}, {0.002, 4.0E-4}, {0, 4.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.002, -2.0E-4}, {0.002, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.002, 0.002}, {0.002, 4.0E-4}, {0, 4.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.002, -2.0E-4}, {0.002, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.002, 0.002}, {0.002, 2.0E-4}, {0, 2.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.002, -2.0E-4}, {0.002, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.0018, 0.002}, {0.0018, 2.0E-4}, {0, 2.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.0018, -2.0E-4}, {0.0018, -0.002}, {0.003, -0.002}, {0.003, 0.002}, {0.0018, 0.002}, {0.0018, 2.0E-4}, {0, 2.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.0018, -2.0E-4}, {0.0018, -0.002}, {0.003310191910713911, -0.002}, {0.003310191910713911, 0.002}, {0.0018, 0.002}, {0.0018, 2.0E-4}, {0, 2.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new double[][]{{0, 0.003}, {0.004, 0.003}, {0.004, -0.003}, {0, -0.003}, {0, -2.0E-4}, {0.0018, -2.0E-4}, {0.0018, -0.002}, {0.0034, -0.002}, {0.0034, 0.002}, {0.0018, 0.002}, {0.0018, 2.0E-4}, {0, 2.0E-4}, {0, 0.003}});
    model.component("comp1").geom("geom1").run("c1");
    model.component("comp1").geom("geom1").feature("c1").set("pos", new double[]{0.0027, 0});
    model.component("comp1").geom("geom1").run("c1");

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").propertyGroup().create("DielectricLoss", "Dielectric losses");
    model.component("comp1").material("mat1").propertyGroup().create("MagneticLosses", "Magnetic losses");
    model.component("comp1").material("mat1").label("Ferrite");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonBis", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonPrim", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("eta_epsilon", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonBis", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonPrim", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("eta_epsilon", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonBis", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("epsilonPrim", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss").set("eta_epsilon", "");
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss")
         .set("epsilonBis", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss")
         .set("epsilonPrim", new String[]{"100000", "0", "0", "0", "100000", "0", "0", "0", "100000"});
    model.component("comp1").material("mat1").propertyGroup("DielectricLoss")
         .set("eta_epsilon", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murPrim", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murBis", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murPrim", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murBis", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murPrim", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murBis", "");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murPrim", "3000");
    model.component("comp1").material("mat1").propertyGroup("MagneticLosses").set("murBis", "0");
    model.component("comp1").material("mat1").set("family", "plastic");

    model.component("comp1").geom("geom1").run();

    model.component("comp1").material("mat1").selection().set(1);
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat2").propertyGroup().create("Enu", "Young's modulus and Poisson's ratio");
    model.component("comp1").material("mat2").propertyGroup().create("linzRes", "Linearized resistivity");
    model.component("comp1").material("mat2").label("Copper");
    model.component("comp1").material("mat2").set("family", "copper");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("relpermeability", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("electricconductivity", new String[]{"5.998e7[S/m]", "0", "0", "0", "5.998e7[S/m]", "0", "0", "0", "5.998e7[S/m]"});
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"17e-6[1/K]", "0", "0", "0", "17e-6[1/K]", "0", "0", "0", "17e-6[1/K]"});
    model.component("comp1").material("mat2").propertyGroup("def").set("heatcapacity", "385[J/(kg*K)]");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat2").propertyGroup("def").set("density", "8960[kg/m^3]");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("thermalconductivity", new String[]{"400[W/(m*K)]", "0", "0", "0", "400[W/(m*K)]", "0", "0", "0", "400[W/(m*K)]"});
    model.component("comp1").material("mat2").propertyGroup("Enu").set("E", "");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("nu", "");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("E", "");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("nu", "");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("E", "110e9[Pa]");
    model.component("comp1").material("mat2").propertyGroup("Enu").set("nu", "0.35");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("rho0", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("alpha", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("Tref", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("rho0", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("alpha", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("Tref", "");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("rho0", "1.72e-8[ohm*m]");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("alpha", "0.0039[1/K]");
    model.component("comp1").material("mat2").propertyGroup("linzRes").set("Tref", "298[K]");
    model.component("comp1").material("mat2").propertyGroup("linzRes").addInput("temperature");
    model.component("comp1").material("mat2").set("family", "copper");
    model.component("comp1").material("mat2").selection().set(2);
    model.component("comp1").material().create("mat3", "Common");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("eta", "Piecewise");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("rho", "Analytic");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("k", "Piecewise");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("cs", "Analytic");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("an1", "Analytic");
    model.component("comp1").material("mat3").propertyGroup("def").func().create("an2", "Analytic");
    model.component("comp1").material("mat3").propertyGroup().create("RefractiveIndex", "Refractive index");
    model.component("comp1").material("mat3").propertyGroup().create("NonlinearModel", "Nonlinear model");
    model.component("comp1").material("mat3").propertyGroup().create("idealGas", "Ideal gas");
    model.component("comp1").material("mat3").propertyGroup("idealGas").func().create("Cp", "Piecewise");
    model.component("comp1").material("mat3").label("Air");
    model.component("comp1").material("mat3").set("family", "air");
    model.component("comp1").material("mat3").propertyGroup("def").func("eta").set("arg", "T");
    model.component("comp1").material("mat3").propertyGroup("def").func("eta")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("eta").set("argunit", "K");
    model.component("comp1").material("mat3").propertyGroup("def").func("eta").set("fununit", "Pa*s");
    model.component("comp1").material("mat3").propertyGroup("def").func("Cp").set("arg", "T");
    model.component("comp1").material("mat3").propertyGroup("def").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat3").propertyGroup("def").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat3").propertyGroup("def").func("rho")
         .set("expr", "pA*0.02897/R_const[K*mol/J]/T");
    model.component("comp1").material("mat3").propertyGroup("def").func("rho").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat3").propertyGroup("def").func("rho").set("fununit", "kg/m^3");
    model.component("comp1").material("mat3").propertyGroup("def").func("rho")
         .set("argunit", new String[]{"Pa", "K"});
    model.component("comp1").material("mat3").propertyGroup("def").func("rho")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "293.15"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("k").set("arg", "T");
    model.component("comp1").material("mat3").propertyGroup("def").func("k")
         .set("pieces", new String[][]{{"200.0", "1600.0", "-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("k").set("argunit", "K");
    model.component("comp1").material("mat3").propertyGroup("def").func("k").set("fununit", "W/(m*K)");
    model.component("comp1").material("mat3").propertyGroup("def").func("cs")
         .set("expr", "sqrt(1.4*R_const[K*mol/J]/0.02897*T)");
    model.component("comp1").material("mat3").propertyGroup("def").func("cs").set("args", new String[]{"T"});
    model.component("comp1").material("mat3").propertyGroup("def").func("cs").set("fununit", "m/s");
    model.component("comp1").material("mat3").propertyGroup("def").func("cs").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat3").propertyGroup("def").func("cs")
         .set("plotargs", new String[][]{{"T", "273.15", "373.15"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("an1").set("funcname", "alpha_p");
    model.component("comp1").material("mat3").propertyGroup("def").func("an1")
         .set("expr", "-1/rho(pA,T)*d(rho(pA,T),T)");
    model.component("comp1").material("mat3").propertyGroup("def").func("an1").set("args", new String[]{"pA", "T"});
    model.component("comp1").material("mat3").propertyGroup("def").func("an1").set("fununit", "1/K");
    model.component("comp1").material("mat3").propertyGroup("def").func("an1")
         .set("argunit", new String[]{"Pa", "K"});
    model.component("comp1").material("mat3").propertyGroup("def").func("an1")
         .set("plotargs", new String[][]{{"pA", "101325", "101325"}, {"T", "273.15", "373.15"}});
    model.component("comp1").material("mat3").propertyGroup("def").func("an2").set("funcname", "muB");
    model.component("comp1").material("mat3").propertyGroup("def").func("an2").set("expr", "0.6*eta(T)");
    model.component("comp1").material("mat3").propertyGroup("def").func("an2").set("args", new String[]{"T"});
    model.component("comp1").material("mat3").propertyGroup("def").func("an2").set("fununit", "Pa*s");
    model.component("comp1").material("mat3").propertyGroup("def").func("an2").set("argunit", new String[]{"K"});
    model.component("comp1").material("mat3").propertyGroup("def").func("an2")
         .set("plotargs", new String[][]{{"T", "200", "1600"}});
    model.component("comp1").material("mat3").propertyGroup("def").set("thermalexpansioncoefficient", "");
    model.component("comp1").material("mat3").propertyGroup("def").set("molarmass", "");
    model.component("comp1").material("mat3").propertyGroup("def").set("bulkviscosity", "");
    model.component("comp1").material("mat3").propertyGroup("def")
         .set("thermalexpansioncoefficient", new String[]{"alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)", "0", "0", "0", "alpha_p(pA,T)"});
    model.component("comp1").material("mat3").propertyGroup("def").set("molarmass", "0.02897[kg/mol]");
    model.component("comp1").material("mat3").propertyGroup("def").set("bulkviscosity", "muB(T)");
    model.component("comp1").material("mat3").propertyGroup("def")
         .set("relpermeability", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat3").propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat3").propertyGroup("def").set("dynamicviscosity", "eta(T)");
    model.component("comp1").material("mat3").propertyGroup("def").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat3").propertyGroup("def")
         .set("electricconductivity", new String[]{"0[S/m]", "0", "0", "0", "0[S/m]", "0", "0", "0", "0[S/m]"});
    model.component("comp1").material("mat3").propertyGroup("def").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat3").propertyGroup("def").set("density", "rho(pA,T)");
    model.component("comp1").material("mat3").propertyGroup("def")
         .set("thermalconductivity", new String[]{"k(T)", "0", "0", "0", "k(T)", "0", "0", "0", "k(T)"});
    model.component("comp1").material("mat3").propertyGroup("def").set("soundspeed", "cs(T)");
    model.component("comp1").material("mat3").propertyGroup("def").addInput("temperature");
    model.component("comp1").material("mat3").propertyGroup("def").addInput("pressure");
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex").set("n", "");
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex").set("ki", "");
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex")
         .set("n", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat3").propertyGroup("RefractiveIndex")
         .set("ki", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat3").propertyGroup("NonlinearModel").set("BA", "(def.gamma+1)/2");
    model.component("comp1").material("mat3").propertyGroup("idealGas").func("Cp").label("Piecewise 2");
    model.component("comp1").material("mat3").propertyGroup("idealGas").func("Cp").set("arg", "T");
    model.component("comp1").material("mat3").propertyGroup("idealGas").func("Cp")
         .set("pieces", new String[][]{{"200.0", "1600.0", "1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4"}});
    model.component("comp1").material("mat3").propertyGroup("idealGas").func("Cp").set("argunit", "K");
    model.component("comp1").material("mat3").propertyGroup("idealGas").func("Cp").set("fununit", "J/(kg*K)");
    model.component("comp1").material("mat3").propertyGroup("idealGas").set("Rs", "R_const/Mn");
    model.component("comp1").material("mat3").propertyGroup("idealGas").set("heatcapacity", "Cp(T)");
    model.component("comp1").material("mat3").propertyGroup("idealGas").set("ratioofspecificheat", "1.4");
    model.component("comp1").material("mat3").propertyGroup("idealGas").set("molarmass", "0.02897");
    model.component("comp1").material("mat3").propertyGroup("idealGas").addInput("temperature");
    model.component("comp1").material("mat3").propertyGroup("idealGas").addInput("pressure");
    model.component("comp1").material("mat3").materialType("nonSolid");
    model.component("comp1").material("mat3").set("family", "air");
    model.component("comp1").material().remove("mat3");
    model.component("comp1").material("mat1").selection().set(1);
    model.component("comp1").material("mat1").set("color", "blue");
    model.component("comp1").material("mat1").propertyGroup("def").set("electricconductivity", new String[]{"0.3"});

    model.study("std1").feature("freq").set("punit", "MHz");
    model.study("std1").feature("freq").set("plist", 1);

    model.component("comp1").physics("mef").create("coil1", "Coil", 2);
    model.component("comp1").physics("mef").feature("coil1").selection().set(2);
    model.component("comp1").physics("mef").feature("coil1").setIndex("materialType", "solid", 0);
    model.component("comp1").physics("mef").feature("coil1").set("N", 1);

    model.sol().create("sol1");
    model.sol("sol1").study("std1");

    model.study("std1").feature("freq").set("notlistsolnum", 1);
    model.study("std1").feature("freq").set("notsolnum", "auto");
    model.study("std1").feature("freq").set("listsolnum", 1);
    model.study("std1").feature("freq").set("solnum", "auto");

    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").feature("st1").set("study", "std1");
    model.sol("sol1").feature("st1").set("studystep", "freq");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").feature("v1").set("control", "freq");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("p1", "Parametric");
    model.sol("sol1").feature("s1").feature().remove("pDef");
    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"freq"});
    model.sol("sol1").feature("s1").feature("p1").set("plistarr", new String[]{"1"});
    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{"MHz"});
    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");
    model.sol("sol1").feature("s1").feature("p1").set("preusesol", "auto");
    model.sol("sol1").feature("s1").feature("p1").set("pdistrib", "off");
    model.sol("sol1").feature("s1").feature("p1").set("plot", "off");
    model.sol("sol1").feature("s1").feature("p1").set("plotgroup", "Default");
    model.sol("sol1").feature("s1").feature("p1").set("probesel", "all");
    model.sol("sol1").feature("s1").feature("p1").set("probes", new String[]{});
    model.sol("sol1").feature("s1").feature("p1").set("control", "freq");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "dDef");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol("sol1").attach("std1");

    model.result().create("pg1", "PlotGroup2D");
    model.result("pg1").label("Magnetic Flux Density Norm (mef)");
    model.result("pg1").set("dataisaxisym", "off");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").set("showlegendsmaxmin", true);
    model.result("pg1").set("data", "dset1");
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("solutionparams", "parent");
    model.result("pg1").feature("surf1").set("expr", "mef.normB");
    model.result("pg1").feature("surf1").set("colortable", "Prism");
    model.result("pg1").feature("surf1").set("colortabletrans", "nonlinear");
    model.result("pg1").feature("surf1").set("colorcalibration", -0.8);
    model.result("pg1").feature("surf1").set("showsolutionparams", "on");
    model.result("pg1").feature("surf1").set("data", "parent");
    model.result("pg1").feature().create("str1", "Streamline");
    model.result("pg1").feature("str1").set("showsolutionparams", "on");
    model.result("pg1").feature("str1").set("solutionparams", "parent");
    model.result("pg1").feature("str1").set("titletype", "none");
    model.result("pg1").feature("str1").set("posmethod", "uniform");
    model.result("pg1").feature("str1").set("udist", 0.03);
    model.result("pg1").feature("str1").set("maxlen", 0.4);
    model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg1").feature("str1").set("inheritcolor", false);
    model.result("pg1").feature("str1").set("showsolutionparams", "on");
    model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg1").feature("str1").set("showsolutionparams", "on");
    model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg1").feature("str1").set("showsolutionparams", "on");
    model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg1").feature("str1").set("showsolutionparams", "on");
    model.result("pg1").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg1").feature("str1").set("data", "parent");
    model.result("pg1").feature("str1").selection().geom("geom1", 1);
    model.result("pg1").feature("str1").selection().set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    model.result("pg1").feature("str1").set("inheritplot", "surf1");
    model.result("pg1").feature("str1").feature().create("col1", "Color");
    model.result("pg1").feature("str1").feature("col1").set("expr", "mef.normB");
    model.result("pg1").feature("str1").feature("col1").set("colortable", "PrismDark");
    model.result("pg1").feature("str1").feature("col1").set("colorlegend", false);
    model.result("pg1").feature("str1").feature("col1").set("colortabletrans", "nonlinear");
    model.result("pg1").feature("str1").feature("col1").set("colorcalibration", -0.8);
    model.result("pg1").feature("str1").feature().create("filt1", "Filter");
    model.result("pg1").feature("str1").feature("filt1").set("expr", "!isScalingSystemDomain");
    model.result("pg1").feature().create("con1", "Contour");
    model.result("pg1").feature("con1").set("showsolutionparams", "on");
    model.result("pg1").feature("con1").set("solutionparams", "parent");
    model.result("pg1").feature("con1").set("expr", "r*mef.Aphi");
    model.result("pg1").feature("con1").set("titletype", "none");
    model.result("pg1").feature("con1").set("number", 10);
    model.result("pg1").feature("con1").set("levelrounding", false);
    model.result("pg1").feature("con1").set("coloring", "uniform");
    model.result("pg1").feature("con1").set("colorlegend", false);
    model.result("pg1").feature("con1").set("color", "custom");
    model.result("pg1").feature("con1")
         .set("customcolor", new double[]{0.3764705955982208, 0.3764705955982208, 0.3764705955982208});
    model.result("pg1").feature("con1").set("resolution", "fine");
    model.result("pg1").feature("con1").set("inheritcolor", false);
    model.result("pg1").feature("con1").set("showsolutionparams", "on");
    model.result("pg1").feature("con1").set("data", "parent");
    model.result("pg1").feature("con1").set("inheritplot", "surf1");
    model.result("pg1").feature("con1").feature().create("filt1", "Filter");
    model.result("pg1").feature("con1").feature("filt1").set("expr", "!isScalingSystemDomain");
    model.result().dataset().create("rev1", "Revolve2D");
    model.result().dataset("rev1").set("startangle", -90);
    model.result().dataset("rev1").set("revangle", 225);
    model.result().dataset("rev1").set("data", "dset1");
    model.result().create("pg2", "PlotGroup3D");
    model.result("pg2").label("Magnetic Flux Density Norm, Revolved Geometry (mef)");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").set("showlegendsmaxmin", true);
    model.result("pg2").set("data", "rev1");
    model.result("pg2").feature().create("vol1", "Volume");
    model.result("pg2").feature("vol1").set("showsolutionparams", "on");
    model.result("pg2").feature("vol1").set("solutionparams", "parent");
    model.result("pg2").feature("vol1").set("expr", "mef.normB");
    model.result("pg2").feature("vol1").set("colortable", "Prism");
    model.result("pg2").feature("vol1").set("colortabletrans", "nonlinear");
    model.result("pg2").feature("vol1").set("colorcalibration", -0.8);
    model.result("pg2").feature("vol1").set("showsolutionparams", "on");
    model.result("pg2").feature("vol1").set("data", "parent");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("solutionparams", "parent");
    model.result("pg2").feature("con1").set("expr", "r*mef.Aphi");
    model.result("pg2").feature("con1").set("titletype", "none");
    model.result("pg2").feature("con1").set("number", 10);
    model.result("pg2").feature("con1").set("levelrounding", false);
    model.result("pg2").feature("con1").set("coloring", "uniform");
    model.result("pg2").feature("con1").set("colorlegend", false);
    model.result("pg2").feature("con1").set("color", "custom");
    model.result("pg2").feature("con1")
         .set("customcolor", new double[]{0.3764705955982208, 0.3764705955982208, 0.3764705955982208});
    model.result("pg2").feature("con1").set("resolution", "fine");
    model.result("pg2").feature("con1").set("inheritcolor", false);
    model.result("pg2").feature("con1").set("showsolutionparams", "on");
    model.result("pg2").feature("con1").set("data", "parent");
    model.result("pg2").feature("con1").set("inheritplot", "vol1");
    model.result("pg2").feature("con1").feature().create("filt1", "Filter");
    model.result("pg2").feature("con1").feature("filt1").set("expr", "!isScalingSystemDomain");
    model.result("pg2").feature("con1").feature("filt1").set("shownodespec", "on");
    model.result().create("pg3", "PlotGroup2D");
    model.result("pg3").label("Electric Potential (mef)");
    model.result("pg3").set("dataisaxisym", "off");
    model.result("pg3").set("frametype", "spatial");
    model.result("pg3").set("showlegendsmaxmin", true);
    model.result("pg3").set("data", "dset1");
    model.result("pg3").feature().create("surf1", "Surface");
    model.result("pg3").feature("surf1").set("showsolutionparams", "on");
    model.result("pg3").feature("surf1").set("solutionparams", "parent");
    model.result("pg3").feature("surf1").set("colortable", "Dipole");
    model.result("pg3").feature("surf1").set("showsolutionparams", "on");
    model.result("pg3").feature("surf1").set("data", "parent");
    model.result("pg3").feature().create("str1", "Streamline");
    model.result("pg3").feature("str1").set("showsolutionparams", "on");
    model.result("pg3").feature("str1").set("solutionparams", "parent");
    model.result("pg3").feature("str1").set("expr", new String[]{"mef.Er", "mef.Ez"});
    model.result("pg3").feature("str1").set("titletype", "none");
    model.result("pg3").feature("str1").set("posmethod", "uniform");
    model.result("pg3").feature("str1").set("udist", 0.02);
    model.result("pg3").feature("str1").set("maxlen", 0.4);
    model.result("pg3").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg3").feature("str1").set("inheritcolor", false);
    model.result("pg3").feature("str1").set("showsolutionparams", "on");
    model.result("pg3").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg3").feature("str1").set("showsolutionparams", "on");
    model.result("pg3").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg3").feature("str1").set("showsolutionparams", "on");
    model.result("pg3").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg3").feature("str1").set("showsolutionparams", "on");
    model.result("pg3").feature("str1").set("maxtime", Double.POSITIVE_INFINITY);
    model.result("pg3").feature("str1").set("data", "parent");
    model.result("pg3").feature("str1").selection().geom("geom1", 1);
    model.result("pg3").feature("str1").selection().set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    model.result("pg3").feature("str1").set("inheritplot", "surf1");
    model.result("pg3").feature("str1").feature().create("col1", "Color");
    model.result("pg3").feature("str1").feature("col1").set("colortable", "DipoleDark");
    model.result("pg3").feature("str1").feature("col1").set("colorlegend", false);

    return model;
  }


  public static void main(String[] args) {
    Model model = run();
  }

}
