/** @file readParameters.h
	@brief reads the inputs from the input files.
*/
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
void readParams();
void readParams() {
    FILE * pFile = fopen ("../inputs/unifiedInput.json" , "r");
    char buffer[65536];
    rapidjson::FileReadStream is(pFile, buffer, sizeof(buffer));
    rapidjson::Document document;
    document.ParseStream<0, rapidjson::UTF8<>, rapidjson::FileReadStream>(is);

    Pr=document["physicalProperties"]["pressure"].GetDouble();
    Temp0=document["initialConditions"]["temperature"].GetDouble();
    if (document["kinetics"]["kineticModel"].GetString()==std::string("Baser")) {
        kinMod=1;
    } else if (document["kinetics"]["kineticModel"].GetString()==std::string("BaserRx")) {
        kinMod=2;
    } else if (document["kinetics"]["kineticModel"].GetString()==std::string("RF-1")) {
        kinMod=3;
    } else {
        std::cout << "kinetic model unknown in QmomKinetics" << std::endl;
        exit(0);
    }
    X_gel=document["kinetics"]["gelPoint"].GetDouble();
    if (kinMod==1 || kinMod==2) {
        A_OH=document["kinetics"]["gellingReaction"]["frequentialFactor"].GetDouble();
        E_OH=document["kinetics"]["gellingReaction"]["activationEnergy"].GetDouble();
        A_W=document["kinetics"]["blowingReaction"]["frequentialFactor"].GetDouble();
        E_W=document["kinetics"]["blowingReaction"]["activationEnergy"].GetDouble();
        OH_0=document["initialConditions"]["concentrations"]["polyol"].GetDouble();
        NCO_0=document["initialConditions"]["concentrations"]["isocyanate"].GetDouble();
        if (document["kinetics"]["useDilution"].GetBool()){
            dilution=1;
        } else {
            dilution=0;
        }
    } else if (kinMod==3) {
        catalyst=document["initialConditions"]["concentrations"]["catalyst"].GetDouble();
    	polyol1_ini=document["initialConditions"]["concentrations"]["polyol1"].GetDouble();
    	polyol2_ini=document["initialConditions"]["concentrations"]["polyol2"].GetDouble();
    	amine_ini=document["initialConditions"]["concentrations"]["amine"].GetDouble();
    	isocyanate1_ini=document["initialConditions"]["concentrations"]["isocyanate1"].GetDouble();
    	isocyanate2_ini=document["initialConditions"]["concentrations"]["isocyanate2"].GetDouble();
    	isocyanate3_ini=document["initialConditions"]["concentrations"]["isocyanate3"].GetDouble();
        OH_0=polyol1_ini+polyol2_ini;
    }
    W_0=document["initialConditions"]["concentrations"]["water"].GetDouble();
    rhoBL=document["physicalProperties"]["blowingAgents"]["PBL"]["density"].GetDouble();
    DH_OH=document["kinetics"]["gellingReaction"]["reactionEnthalpy"].GetDouble();
    DH_W=document["kinetics"]["blowingReaction"]["reactionEnthalpy"].GetDouble();
    C_Poly=document["physicalProperties"]["polymer"]["heatCapacity"].GetDouble();
    C_CO2=document["physicalProperties"]["blowingAgents"]["CO2"]["heatCapacityInLiquidPhase"].GetDouble();
    C_BG=document["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInGaseousPhase"].GetDouble();
    C_BL=document["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInLiquidPhase"].GetDouble();
    lambda=document["physicalProperties"]["blowingAgents"]["PBL"]["evaporationHeat"].GetDouble();
    if (document["physicalBlowingAgent"].GetString()==std::string("n-pentane")) {
        phBL=1;
    } else if (document["physicalBlowingAgent"].GetString()==std::string("R11")) {
        phBL=2;
    } else {
        std::cout << "unknown blowing agent (solubility)" << std::endl;
        exit(0);
    }
    if (document["physicalProperties"]["polymer"]["polymerDensityModel"].GetString()==std::string("nanotools")) {
        denMod=1;
    } else if (document["physicalProperties"]["polymer"]["polymerDensityModel"].GetString()==std::string("constant")) {
        denMod=2;
        rhoPoly=document["physicalProperties"]["polymer"]["density"].GetDouble();
    } else if (document["physicalProperties"]["polymer"]["polymerDensityModel"].GetString()==std::string("pcsaft")) {
        denMod=3;
    }
    M_CO2=document["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"].GetDouble()*1e3;
    M_B=document["physicalProperties"]["blowingAgents"]["PBL"]["molarMass"].GetDouble()*1e3;
    M_NCO=document["physicalProperties"]["polymer"]["molarMassNCO"].GetDouble()*1e3;
    M_air=document["physicalProperties"]["air"]["molarMass"].GetDouble()*1e3;
    double cb=document["initialConditions"]["concentrations"]["blowingAgents"]["PBL"].GetDouble();
    L0=cb*M_B/rhoPoly/1000.0;
    double cc=document["initialConditions"]["concentrations"]["blowingAgents"]["CO2"].GetDouble();
    CO2_0=cc*M_CO2/rhoPoly/1000.0;
    surfaceTension=document["physicalProperties"]["surfaceTension"].GetDouble();
    sig=document["initialConditions"]["bubbleRadiusDeviation"].GetDouble();
    init_size=document["initialConditions"]["bubbleRadius"].GetDouble()*2;
    NN=document["initialConditions"]["numberBubbleDensity"].GetDouble();
    abs_err=document["QmomKinetics"]["absoluteTolerance"].GetDouble();
    rel_err=document["QmomKinetics"]["relativeTolerance"].GetDouble();
    dt=document["QmomKinetics"]["timeStep"].GetDouble();
    tend=document["QmomKinetics"]["endTime"].GetDouble();
    bubbleMode=document["QmomKinetics"]["bubbleMode"].GetString();
    if (W_0<1e-8) {
        W_0=-1.0;
    }
    apparentViscosity=document["physicalProperties"]["ModenaFoamViscosityModel"].GetBool();
}
