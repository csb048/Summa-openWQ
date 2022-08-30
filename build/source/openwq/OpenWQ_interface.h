/**
 * This is the C interface for SUMMA, these are the functions that are called 
 * by SUMMA and the iso bindings. 
 * These are only their definition and their actual implementation is in
 * OpenWQ_hydrolink.cpp 
 */

#ifdef __cplusplus
extern "C" { 
    class ClassWQ_OpenWQ;
    typedef ClassWQ_OpenWQ CLASSWQ_OPENWQ;
    #else
    typedef struct CLASSWQ_OPENWQ CLASSWQ_OPENWQ;
    #endif

    // Create OpenWQ Object
    CLASSWQ_OPENWQ* create_openwq();

    // Delete OpenWQ Object
    void delete_openwq(CLASSWQ_OPENWQ* openWQ);

    // OpenWQ initalization method
    int openwq_decl(CLASSWQ_OPENWQ *openWQ,int numHRU, int num_layers_canopy, int num_layers_matricHead, int num_layers_aquifer, int num_layers_volFracWat, int y_direction);

    int openwq_run_time_start(CLASSWQ_OPENWQ *openWQ, int numHRU, int simtime_summa[],
        double soilMoisture[], double soilTemp[], double airTemp[], double SWE_vol[], double canopyWat[], double matricHead_vol[], double aquiferStorage[]);

    // OpenWQ run functions, this function decides which C++ code to call
    int openwq_run_space(CLASSWQ_OPENWQ *openWQ, int simtime_summa[], int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r, double wmass_source);

    int openwq_run_space_in(CLASSWQ_OPENWQ *openWQ, int simtime_summa[], int recipient, int ix_r, int iy_r, int iz_r, double wflux_s2r);

    int openwq_run_time_end(CLASSWQ_OPENWQ *openWQ, int simtime_summa[]);

    #ifdef __cplusplus
}
#endif