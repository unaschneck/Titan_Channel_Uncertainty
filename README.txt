Titan river analysis pipeline (originally written by Frank Verdi 2024 and refined by Una Schneck 2026)

Sinuosity2 -- 

there was an error in the ordering of the dimensionless coordinates (yd = file.geometry.center_dimensionless(:,1);xd = file.geometry.center_dimensionless(:,2); should be swapped

excluded endpoints in original code 

There appears to be some SAR data that is being identified but not the right dimension size to be used for extracting brightness and is not 
called out in the table in the draft. Why is this?