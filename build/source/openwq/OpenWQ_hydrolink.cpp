// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"


// Constructor
// initalize numHRUs value
ClassWQ_OpenWQ::ClassWQ_OpenWQ() {}

// Deconstructor
ClassWQ_OpenWQ::~ClassWQ_OpenWQ() {}

time_t ClassWQ_OpenWQ::convert_time(
    int year, 
    int month, 
    int day, 
    int hour, 
    int minute) {

    std::time_t sim_time;
    std::tm tm{};
    tm.tm_year = year - 1900; // -1900 is needed to get the conversion to produce the correct output
    tm.tm_mon = month - 1;
    tm.tm_hour = hour;
    tm.tm_mday = day;
    tm.tm_min = minute;
    sim_time = timegm(&tm);

    return sim_time;
}

int ClassWQ_OpenWQ::decl(
    int num_HRU,                // num HRU
    int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
    int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
    int nSoil_2openwq,        // num layers of snoil (variable)
    int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    int nYdirec_2openwq){           // num of layers in y-dir (set to 1 because not used in summa)
    
    OpenWQ_hostModelconfig_ref = new OpenWQ_hostModelconfig();
    OpenWQ_couplercalls_ref = new OpenWQ_couplercalls();
    OpenWQ_json_ref = new OpenWQ_json();
    OpenWQ_wqconfig_ref = new OpenWQ_wqconfig();
    OpenWQ_units_ref = new OpenWQ_units();
    OpenWQ_utils_ref = new OpenWQ_utils();
    OpenWQ_readjson_ref = new OpenWQ_readjson();
    OpenWQ_initiate_ref = new OpenWQ_initiate();
    OpenWQ_watertransp_ref = new OpenWQ_watertransp();
    OpenWQ_chem_ref = new OpenWQ_chem();
    OpenWQ_extwatflux_ss_ref = new OpenWQ_extwatflux_ss();
    OpenWQ_output_ref = new OpenWQ_output();
    
    this->num_HRU = num_HRU;

    if (OpenWQ_hostModelconfig_ref->HydroComp.size()==0) {

        // Compartment names
        // Make sure to use capital letters for compartment names
        OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"SCALARCANOPYWAT",num_HRU,nYdirec_2openwq,nCanopy_2openwq));      // Canopy
        OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(1,"ILAYERVOLFRACWAT",num_HRU,nYdirec_2openwq,nSnow_2openwq + nSoil_2openwq)); // Soil + Snow
        OpenWQ_hostModelconfig_ref->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(2,"SCALARAQUIFER",num_HRU,nYdirec_2openwq,nAquifer_2openwq));       // GW
        

        OpenWQ_vars_ref = new OpenWQ_vars(OpenWQ_hostModelconfig_ref->HydroComp.size());

        // External fluxes
        // Make sure to use capital letters for external fluxes
        OpenWQ_hostModelconfig_ref->HydroExtFlux.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"PRECIP",num_HRU,nYdirec_2openwq,1));

        // Dependencies
        // to expand BGC modelling options
        OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"SM",num_HRU,nYdirec_2openwq,1));
        OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(1,"Tair_K",num_HRU,nYdirec_2openwq,1));
        OpenWQ_hostModelconfig_ref->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(2,"Tsoil_K",num_HRU,nYdirec_2openwq,1));

        // Master Json
        OpenWQ_wqconfig_ref->OpenWQ_masterjson = "openWQ_master.json";


        OpenWQ_couplercalls_ref->InitialConfig(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,                // create OpenWQ_json object
            *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
            *OpenWQ_units_ref,               // functions for unit conversion
            *OpenWQ_utils_ref,                // utility methods/functions
            *OpenWQ_readjson_ref,            // read json files
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,            // initiate modules
            *OpenWQ_watertransp_ref,         // transport modules
            *OpenWQ_chem_ref,                // biochemistry modules
            *OpenWQ_extwatflux_ss_ref,       // sink and source modules)
            *OpenWQ_output_ref);
            
    }
    return 0;
}

// soilMoist_depVar does not have a value - it is passed as 0
int ClassWQ_OpenWQ::run_time_start(
    int numHRU, 
    int nSnow_2openwq, 
    int nSoil_2openwq,
    int simtime_summa[], 
    double soilMoist_depVar[], 
    double soilTemp_K_depVar[], 
    double airTemp_K_depVar[],
    double sweWatVol_stateVar[], 
    double canopyWat[], 
    double soilWatVol_stateVar[], 
    double aquiferStorage[]) {

    time_t simtime = convert_time(simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4]);

    for (int i = 0; i < numHRU; i++) {
        // Updating Chemistry dependencies
        (*OpenWQ_hostModelconfig_ref->dependVar)[0](i,0,0) = soilMoist_depVar[i]; 
        (*OpenWQ_hostModelconfig_ref->dependVar)[1](i,0,0) = airTemp_K_depVar[i];
        (*OpenWQ_hostModelconfig_ref->dependVar)[2](i,0,0) = soilTemp_K_depVar[i];
        // Updating water volumes
        //(*OpenWQ_hostModelconfig_ref->waterVol_hydromodel)[0](i,0,0) = sweWatVol_stateVar[i];
        (*OpenWQ_hostModelconfig_ref->waterVol_hydromodel)[0](i,0,0) = canopyWat[i];
        (*OpenWQ_hostModelconfig_ref->waterVol_hydromodel)[1](i,0,0) = soilWatVol_stateVar[i];
        (*OpenWQ_hostModelconfig_ref->waterVol_hydromodel)[2](i,0,0) = aquiferStorage[i];
    }

    // *OpenWQ_hostModelconfig_ref.time_step = 5;

    OpenWQ_couplercalls_ref->RunTimeLoopStart(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
        *OpenWQ_units_ref,               // functions for unit conversion
        *OpenWQ_utils_ref,                // utility methods/functions
        *OpenWQ_readjson_ref,            // read json files
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,            // initiate modules
        *OpenWQ_watertransp_ref,         // transport modules
        *OpenWQ_chem_ref,                // biochemistry modules
        *OpenWQ_extwatflux_ss_ref,          // sink and source modules)
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime);

    return 0;
}

int ClassWQ_OpenWQ::run_space(
    int simtime_summa[], 
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r, double wmass_source) {

    // Convert Fortran Index to C++ index
    ix_s -= 1; iy_s -= 1; iz_s -= 1;
    ix_r -= 1; iy_r -= 1; iz_r -= 1;

   
    time_t simtime = convert_time(
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4]);
    
    OpenWQ_couplercalls_ref->RunSpaceStep(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
        *OpenWQ_units_ref,               // functions for unit conversion
        *OpenWQ_utils_ref,                // utility methods/functions
        *OpenWQ_readjson_ref,            // read json files
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,            // initiate modules
        *OpenWQ_watertransp_ref,         // transport modules
        *OpenWQ_chem_ref,                // biochemistry modules
        *OpenWQ_extwatflux_ss_ref,       // sink and source modules
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);

    return 0;
}

int ClassWQ_OpenWQ::run_space_in(
    int simtime_summa[],
    std::string source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r) {

    // Convert Fortran Index to C++ index
    ix_r -= 1; iy_r -= 1; iz_r -= 1;
    
    time_t simtime = convert_time(
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4]);

     OpenWQ_couplercalls_ref->RunSpaceStep_IN(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_watertransp_ref,
        *OpenWQ_chem_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

    return 0;
}

int ClassWQ_OpenWQ::run_time_end(
    int simtime_summa[]) {
    
    time_t simtime = convert_time(
        simtime_summa[0], 
        simtime_summa[1], 
        simtime_summa[2], 
        simtime_summa[3], 
        simtime_summa[4]);


    OpenWQ_couplercalls_ref->RunTimeLoopEnd(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,            // create OpenWQ_wqconfig object
        *OpenWQ_units_ref,               // functions for unit conversion
        *OpenWQ_utils_ref,                // utility methods/functions
        *OpenWQ_readjson_ref,            // read json files
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,            // initiate modules
        *OpenWQ_watertransp_ref,         // transport modules
        *OpenWQ_chem_ref,                // biochemistry modules
        *OpenWQ_extwatflux_ss_ref,          // sink and source modules)
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime);

    return 0;
}

