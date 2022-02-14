echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent\ntbin\win64\tell.exe" DESKTOP-R5L45J1 61287 CLEANUP_EXITING
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 8620) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 22340) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 25800)
del "C:\Users\rylab\Desktop\microfluidics-inertial-lift\microfluidic-simulation_files\user_files\cleanup-fluent-DESKTOP-R5L45J1-22340.bat"
