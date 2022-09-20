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
    int nRunoff_2openwq,      // num layers of runoff (fixed to 1)
    int nAquifer_2openwq,     // num layers of aquifer (fixed to 1)
    int nYdirec_2openwq){            // num of layers in y-dir (set to 1 because not used in summa)

    return openWQ->decl(
        hruCount, 
        nCanopy_2openwq, 
        nSnow_2openwq, 
        nSoil_2openwq, 
        nRunoff_2openwq,
        nAquifer_2openwq, 
        nYdirec_2openwq);

}


int openwq_run_time_start(
    ClassWQ_OpenWQ *openWQ,
    bool last_hru_flag, 
    int hru_index, 
    int nSnow_2openwq, 
    int nSoil_2openwq, 
    int simtime_summa[], 
    double soilMoist_depVar[], 
    double soilTemp_K_depVar[], 
    double airTemp_K_depVar,
    double sweWatVol_stateVar[], 
    double canopyWatVol_stateVar, 
    double soilWatVol_stateVar[], 
    double aquiferWatVol_stateVar) {
    
    return openWQ->run_time_start(
        last_hru_flag,
        hru_index, 
        nSnow_2openwq, 
        nSoil_2openwq, 
        simtime_summa, 
        soilMoist_depVar, 
        soilTemp_K_depVar, 
        airTemp_K_depVar, 
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
    char* source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r, 
    double wflux_s2r) {
    
    // convert source_EWF_name to string
    std::string source_EWF_name_str(source_EWF_name);

    return openWQ->run_space_in(
        simtime_summa,
        source_EWF_name_str,
        recipient, ix_r, iy_r, iz_r, 
        wflux_s2r);
}


int openwq_run_time_end(
    ClassWQ_OpenWQ *openWQ, 
    int simtime_summa[]) {

    return openWQ->run_time_end(
        simtime_summa);
}
