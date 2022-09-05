#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"
/**
 * Below is the implementation of the C interface for SUMMA. When Summa calls a function 
 * the functions below are the ones that are invoked first. 
 * The openWQ object is then passed from Fortran to these functions so that the OpenWQ object
 * can be called. The openWQ object methods are defined above.
 */
// Interface functions to create Object
CLASSWQ_OPENWQ* create_openwq() {
    return new ClassWQ_OpenWQ();
}

void delete_openwq(CLASSWQ_OPENWQ* openWQ) {
    delete openWQ;
}

int openwq_decl(ClassWQ_OpenWQ *openWQ, int numHRU, int num_layers_canopy, int num_layers_matricHead, int num_layers_aquifer, int num_layers_volFracWat, int y_direction) {
    return openWQ->decl(numHRU, num_layers_canopy, num_layers_matricHead, num_layers_aquifer, num_layers_volFracWat, y_direction);
}


int openwq_run_time_start(ClassWQ_OpenWQ *openWQ, int numHRU, int maxNumLayers_snow, int maxNumLayers_soil, 
    int simtime_summa[], double soilMoisture[], double soilTemp[], double airTemp[],
    double SWE_vol[], double canopyWat_vol[], double matricHead_vol[], double aquiferStorage_vol[]) {
    
    return openWQ->run_time_start(numHRU, maxNumLayers_snow, maxNumLayers_soil, simtime_summa, 
        soilMoisture, soilTemp, airTemp, SWE_vol, canopyWat_vol, matricHead_vol, aquiferStorage_vol);
}


int openwq_run_space(ClassWQ_OpenWQ *openWQ, int simtime_summa[], int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r, double wmass_source) {

    return openWQ->run_space(simtime_summa, source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r, wflux_s2r, wmass_source);
}

int openwq_run_space_in(ClassWQ_OpenWQ *openWQ, int simtime_summa[], int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r) {
    
    return openWQ->run_space_in(simtime_summa, recipient, ix_r, iy_r, iz_r, wflux_s2r);
}


int openwq_run_time_end(ClassWQ_OpenWQ *openWQ, int simtime_summa[]) {

    return openWQ->run_time_end(simtime_summa);
}
