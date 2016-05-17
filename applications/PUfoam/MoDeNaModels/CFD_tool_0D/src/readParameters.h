/** @file readParameters.h
	@brief reads the inputs from the input files.
*/
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>>
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
void readParams();
void readParams() {
    FILE * pFile = fopen ("../unifiedInput.json" , "r");
    char buffer[65536];
    rapidjson::FileReadStream is(pFile, buffer, sizeof(buffer));
    rapidjson::Document document;
    document.ParseStream<0, rapidjson::UTF8<>, rapidjson::FileReadStream>(is);

    Pr=document["physicalProperties"]["pressure"].GetDouble();
    Temp0=document["initialConditions"]["temperature"].GetDouble();
    A_OH=document["kinetics"]["gellingReaction"]["frequentialFactor"].GetDouble();
    E_OH=document["kinetics"]["gellingReaction"]["activationEnergy"].GetDouble();
    OH_0=document["initialConditions"]["concentrations"]["polyol"].GetDouble();
    NCO_0=document["initialConditions"]["concentrations"]["isocyanate"].GetDouble();
    W_0=document["initialConditions"]["concentrations"]["water"].GetDouble();
    A_W=document["kinetics"]["blowingReaction"]["frequentialFactor"].GetDouble();
    E_W=document["kinetics"]["blowingReaction"]["activationEnergy"].GetDouble();
    rhoPoly=document["physicalProperties"]["polymer"]["density"].GetDouble();
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
    denMod=2; //TODO implement
    if (document["kinetics"]["kineticModel"].GetString()==std::string("Baser")) {
        kinMod=1;
    } else if (document["kinetics"]["kineticModel"].GetString()==std::string("BaserRx")) {
        kinMod=2;
    } else if (document["kinetics"]["kineticModel"].GetString()==std::string("modena")) {
        kinMod=3;
    } else {
        std::cout << "kinetic model unknown in QmomKinetics" << std::endl;
        exit(0);
    }
    if (document["kinetics"]["useDilution"].GetBool()){
        dilution=1;
    } else {
        dilution=0;
    }
    M_CO2=document["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"].GetDouble()*1e3;
    M_B=document["physicalProperties"]["blowingAgents"]["PBL"]["molarMass"].GetDouble()*1e3;
    M_NCO=document["physicalProperties"]["polymer"]["molarMassNCO"].GetDouble()*1e3;
    M_air=document["physicalProperties"]["air"]["molarMass"].GetDouble()*1e3;
    double H=document["physicalProperties"]["blowingAgents"]["CO2"]["solubility"].GetDouble();
    CO2_D=H*Pr*M_CO2/rhoPoly;
    double cb=document["initialConditions"]["concentrations"]["blowingAgents"]["PBL"].GetDouble();
    L0=cb*M_B/rhoPoly;
    double cc=document["initialConditions"]["concentrations"]["blowingAgents"]["CO2"].GetDouble();
    CO2_0=cc*M_CO2/rhoPoly;
    surfaceTension=document["physicalProperties"]["surfaceTension"].GetDouble();
    sig=document["initialConditions"]["bubbleRadiusDeviation"].GetDouble();
    init_size=document["initialConditions"]["bubbleRadius"].GetDouble()*2;
    NN=document["initialConditions"]["numberBubbleDensity"].GetDouble();
    abs_err=document["QmomKinetics"]["absoluteTolerance"].GetDouble();
    rel_err=document["QmomKinetics"]["relativeTolerance"].GetDouble();
    dt=document["QmomKinetics"]["timeStep"].GetDouble();
    tend=document["QmomKinetics"]["endTime"].GetDouble();
    if (W_0<1e-8) {
        W_0=-1.0;
    }
}
