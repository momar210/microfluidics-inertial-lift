echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent\ntbin\win64\tell.exe" DESKTOP-R5L45J1 58401 CLEANUP_EXITING
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 13640) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 14376) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 3660) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 12624) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 18240) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 6228)
del "C:\Users\rylab\Desktop\microfluidics-inertial-lift\microfluidic-simulation_files\user_files\cleanup-fluent-DESKTOP-R5L45J1-18240.bat"
