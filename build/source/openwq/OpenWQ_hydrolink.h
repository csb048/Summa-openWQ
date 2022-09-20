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

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED

#include "OpenWQ_couplercalls.h"
#include "OpenWQ_global.h"
#include "OpenWQ_readjson.h"
#include "OpenWQ_initiate.h"
#include "OpenWQ_chem.h"
#include "OpenWQ_watertransp.h"
#include "OpenWQ_extwatflux_ss.h"
#include "OpenWQ_units.h"
#include "OpenWQ_utils.h"
#include "OpenWQ_solver.h"
#include "OpenWQ_output.h"
#include <iostream>
#include <time.h>
#include <vector>

// Global Indexes for Compartments
  inline int canopy_index_openwq    = 0;
  inline int snow_index_openwq      = 1;
  inline int runoff_index_openwq    = 2;
  inline int soil_index_openwq      = 3;
  inline int aquifer_index_openwq   = 4;

class ClassWQ_OpenWQ
{

    // Instance Variables
    private:

        OpenWQ_couplercalls *OpenWQ_couplercalls_ref;
        OpenWQ_hostModelconfig *OpenWQ_hostModelconfig_ref;
        OpenWQ_json *OpenWQ_json_ref;
        OpenWQ_wqconfig *OpenWQ_wqconfig_ref;
        OpenWQ_units *OpenWQ_units_ref;
        OpenWQ_utils *OpenWQ_utils_ref;
        OpenWQ_readjson *OpenWQ_readjson_ref;
        OpenWQ_vars *OpenWQ_vars_ref;
        OpenWQ_initiate *OpenWQ_initiate_ref;
        OpenWQ_watertransp *OpenWQ_watertransp_ref;
        OpenWQ_chem *OpenWQ_chem_ref;            
        OpenWQ_extwatflux_ss *OpenWQ_extwatflux_ss_ref;
        OpenWQ_solver *OpenWQ_solver_ref;
        OpenWQ_output *OpenWQ_output_ref;

        int num_HRU;
        const float *hru_area;

    // Constructor
    public:
        ClassWQ_OpenWQ();
        ~ClassWQ_OpenWQ();
    
    // Methods
    void printNum() {
        std::cout << "num = " << this->num_HRU << std::endl;
    }

    int decl(
        int num_HRU,                // num HRU
        int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
        int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
        int nSoil_2openwq,        // num layers of snoil (variable)
        int nRunoff_2openwq,      // num layers of runoff (fixed to 1)
        int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
        int nYdirec_2openwq);           // num of layers in y-dir (set to 1 because not used in summa)

    int run_time_start(
        bool last_hru_flag,
        int hru_index, 
        int nSnow_2openwq, 
        int nSoil_2openwq, 
        int simtime_summa[],
        double soilMoist_depVar[], 
        double soilTemp_K_depVar[], 
        double airTemp_K_depVar,
        double sweWatVol_stateVar[], 
        double canopyWat, 
        double soilWatVol_stateVar[], 
        double aquiferStorage);

    int run_space(
        int simtime_summa[], 
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r, 
        double wflux_s2r, double wmass_source);

    int run_space_in(
        int simtime_summa[],
        std::string source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r, 
        double wflux_s2r);

    int run_time_end(
        int simtime_summa[]);

    time_t convert_time(
        int year, 
        int month, 
        int day, 
        int hour, 
        int minute);

};
#endif