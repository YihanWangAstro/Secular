
./secular arg1 arg2 arg3 arg4

/*
arg1 : input file path
arg2 : start id of task in initial file
arg3 : end id of task in initial file
arg4 : name of output directories
*/

initial file format:

task_id   end_time   output_time_interval  Oct GR GW S_inL_out LL  orbit-parameters(like before).......

note:
1.if output_time_interval = 0, no trajectories will be outputed. 
2.GW is the final semi-major axis to stop the simulation(i.e. GW=0.01 means that once 'a' shrink to 0.01 au, the integration will stop). if GW = 0, the GW will be turned off.   !!!
3.the final state of the integration will always be recored.
4. Don't add '.' after the numbers in the initial files.
