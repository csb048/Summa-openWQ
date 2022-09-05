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

int openwq_decl(
    ClassWQ_OpenWQ *openWQ, 
    int hruCount,              // num HRU
    int nCanopy_2openwq,      // num layers of canopy (fixed to 1)
    int nSnow_2openwq,        // num layers of snow (fixed to max of 5 because it varies)
    int nSoil_2openwq,        // num layers of snoil (variable)
    int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    int nYdirec_2openwq){            // num of layers in y-dir (set to 1 because not used in summa)

    return openWQ->decl(
        hruCount, 
        nCanopy_2openwq, 
        nSnow_2openwq, 
        nSoil_2openwq, 
        nAquifer_2openwq, 
        nYdirec_2openwq);

}


int openwq_run_time_start(
    ClassWQ_OpenWQ *openWQ, 
    int numHRU, 
    int nSnow_2openwq, 
    int nSoil_2openwq, 
    int simtime_summa[], 
    double soilMoist_depVar[], 
    double soilTemp_depVar[], 
    double airTemp_depVar[],
    double sweWatVol_stateVar[], 
    double canopyWatVol_stateVar[], 
    double soilWatVol_stateVar[], 
    double aquiferWatVol_stateVar[]) {
    
    return openWQ->run_time_start(
        numHRU, 
        nSnow_2openwq, 
        nSoil_2openwq, 
        simtime_summa, 
        soilMoist_depVar, 
        soilTemp_depVar, 
        airTemp_depVar, 
        sweWatVol_stateVar, 
        canopyWatVol_stateVar, 
        soilWatVol_stateVar, 
        aquiferWatVol_stateVar);
}


int openwq_run_space(
    ClassWQ_OpenWQ *openWQ, 
    int simtime_summa[], 
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r, double wmass_source) {

    return openWQ->run_space(
        simtime_summa, 
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r, 
        wflux_s2r, wmass_source);
}

int openwq_run_space_in(
    ClassWQ_OpenWQ *openWQ, 
    int simtime_summa[], 
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r) {
    
    return openWQ->run_space_in(
        simtime_summa, 
        recipient, ix_r, iy_r, iz_r, 
        wflux_s2r);
}


int openwq_run_time_end(
    ClassWQ_OpenWQ *openWQ, 
    int simtime_summa[]) {

    return openWQ->run_time_end(
        simtime_summa);
}
