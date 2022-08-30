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

        int numHRU;
        const float *hru_area;

    // Constructor
    public:
        ClassWQ_OpenWQ();
        ~ClassWQ_OpenWQ();
    
    // Methods
    void printNum() {
        std::cout << "num = " << this->numHRU << std::endl;
    }

    int decl(int numHRU, int num_layers_canopy, int num_layers_matricHead, int num_layers_aquifer, int num_layers_volFracWat, int y_direction);

    int run_time_start(int numHRU, int simtime_summa[], 
        double soilMoisture[], double soilTemp[], double airTemp[],
        double SWE_vol[], double canopyWat[], double matricHead_vol[], double aquiferStorage[]);

    int run_space(int simtime_summa[], int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r, double wmass_source);

    int run_space_in(int simtime_summa[], int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r);

    int run_time_end(int simtime_summa[]);

    time_t convert_time(int year, int month, int day, int hour, int minute);

};
#endif