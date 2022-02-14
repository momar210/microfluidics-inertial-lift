/* This file generated automatically. */ 
/* Do not modify. */ 
#include "udf.h" 
#include "prop.h" 
#include "dpm.h" 
extern DEFINE_DPM_BODY_FORCE(inertial_lift, p, i);
extern DEFINE_PROFILE(inlet_x_velocity, thread, position);
extern DEFINE_ON_DEMAND(on_demand_max_velocity_calculation);
__declspec(dllexport) UDF_Data udf_data[] = { 
{"inertial_lift", (void (*)(void))inertial_lift, UDF_TYPE_DPM_BODY_FORCE},
{"inlet_x_velocity", (void (*)(void))inlet_x_velocity, UDF_TYPE_PROFILE},
{"on_demand_max_velocity_calculation", (void (*)(void))on_demand_max_velocity_calculation, UDF_TYPE_ON_DEMAND},
}; 
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); 
#include "version.h" 
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) 
{ 
*major = RampantReleaseMajor; 
*minor = RampantReleaseMinor; 
*revision = RampantReleaseRevision; 
} 
