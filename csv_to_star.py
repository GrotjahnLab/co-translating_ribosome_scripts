import pandas as pd
import numpy as np
import glob


header = '''

data_

loop_
_rlnMagnification #1
_rlnDetectorPixelSize #2
_rlnCoordinateX #3
_rlnCoordinateY #4
_rlnCoordinateZ #5
_rlnAngleRot #6
_rlnAngleTilt #7
_rlnAnglePsi #8
_rlnImageName #9
_rlnCtfImage #10
_rlnRandomSubset #11
_rlnPixelSize #12
_rlnVoltage #13
_rlnSphericalAberration #14
_rlnMicrographName #15
_blue #16
_yellow #17
_red #18
_ribo_exit_tunnel_membrane_distance #19
'''

csvfile_list = glob.glob("*.csv")
print(csvfile_list)
for csvfile in csvfile_list:
    df = pd.read_csv(csvfile)
    df = df[df['ribo_exit_tunnel-OMM_dist'] <= 95]
    starfile = csvfile[0:-4] + '_filtered.star'
    #starfile = csvfile[0:-4] + '.star'
    with open(starfile,'w') as outfile:
        outfile.write(header)
        for i in df.index:
            starlines = "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %s %s %.6f %.6f %.6f %.6f %s %.6f %.6f %.6f %.6f\n"%(
                df.loc[i,'rlnMagnification'], df.loc[i,'rlnDetectorPixelSize'],
                df.loc[i,'rlnCoordinateX'], df.loc[i,'rlnCoordinateY'], df.loc[i,'rlnCoordinateZ'],
                df.loc[i, 'rlnAngleRot'], df.loc[i, 'rlnAngleTilt'], df.loc[i, 'rlnAnglePsi'],
                df.loc[i, 'rlnImageName'], df.loc[i, 'rlnCtfImage'], df.loc[i, 'rlnRandomSubset'],
                df.loc[i, 'rlnPixelSize'], df.loc[i, 'rlnVoltage'], df.loc[i, 'rlnSphericalAberration'], df.loc[i, 'rlnMicrographName'], 
                df.loc[i, 'blue'], df.loc[i, 'yellow'], df.loc[i, 'red'], df.loc[i, 'ribo_exit_tunnel-OMM_dist'])
            outfile.write(starlines)
