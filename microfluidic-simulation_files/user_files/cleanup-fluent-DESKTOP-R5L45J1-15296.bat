echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v212\fluent\ntbin\win64\tell.exe" DESKTOP-R5L45J1 59983 CLEANUP_EXITING
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 27876) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 15296) 
if /i "%LOCALHOST%"=="DESKTOP-R5L45J1" (%KILL_CMD% 10720)
del "C:\Users\rylab\Desktop\microfluidics-inertial-lift\microfluidic-simulation_files\user_files\cleanup-fluent-DESKTOP-R5L45J1-15296.bat"
