import pandas as pd
import json
import numpy as np
import re
import os

def shift_centres(sourceinfofile,calfactorsfile,variablesfile,YSOcomparefile,protocatalogue='config/protocat.txt',diskcatalogue='config/diskcat.txt'):
    region = sourceinfofile.split('/')[-1].split('_')[0]
    with open('config/centre_offsets.json') as offsetsfile:
        offsets = json.load(offsetsfile)
        dx,dy = offsets[region][0],offsets[region][1]

    proto_cat     = pd.read_csv(protocatalogue,names=['index','ra','dec','class','cat'],delim_whitespace=True)
    disk_cat      = pd.read_csv(diskcatalogue,names=['index','ra','dec','class','cat'],delim_whitespace=True)
    disk_cat_D    = disk_cat[disk_cat['class'] == 'D']

    #
    # Deal with the sourceinfo file first - we need to update RA, DEC,proto_dist and disk_dist
    #
   
    if os.path.exists(sourceinfofile): 
        sourceinfo = pd.read_csv(sourceinfofile,delim_whitespace=True)
        # RA should be "-". A positive value indicates a DROP in RA, i.e. to the west (looked at individual imags to make sure this was right)
        sourceinfo['RA']  = sourceinfo['RA']-dx/3600.0
        # DEC should be "+". A positive value indicates an INCREASE in DEC, i.e. to the north
        sourceinfo['DEC'] = sourceinfo['DEC']+dy/3600.0

        # Now figure out the closest protostar distances
        protodists   = []
        protoclasses =[]
        for eachsourceRA,eachsourceDEC in zip(sourceinfo['RA'],sourceinfo['DEC']):
            distances = []
            classes   = []
            for eachProtoRA,eachProtoDEC,eachProtoClass in zip(proto_cat['ra'],proto_cat['dec'],proto_cat['class']):
                average_dec = np.average([eachsourceDEC,eachProtoDEC])
                distances.append(3600*np.sqrt(((eachsourceRA-eachProtoRA)*np.cos(average_dec*np.pi/180.0))**2+(eachsourceDEC-eachsourceDEC)**2))
                classes.append(eachProtoClass)
            mindist = np.array(distances)[np.argmin(np.array(distances))]
            protodists.append(mindist)
            protoclasses.append(np.array(classes)[np.argmin(np.array(distances))])

        diskdists   = []
        diskclasses = []
        for eachsourceRA,eachsourceDEC in zip(sourceinfo['RA'],sourceinfo['DEC']):
            distances = []
            classes   = []
            for eachDiskRA,eachDiskDEC,eachDiskClass in zip(disk_cat['ra'],disk_cat['dec'],disk_cat['class']):
                average_dec = np.average([eachsourceDEC,eachDiskDEC])
                distances.append(3600*np.sqrt(((eachsourceRA-eachDiskRA)*np.cos(average_dec*np.pi/180.0))**2+(eachsourceDEC-eachsourceDEC)**2))
                classes.append(eachDiskClass)
            mindist = np.array(distances)[np.argmin(np.array(distances))]
            diskdists.append(mindist)
            diskclasses.append(np.array(classes)[np.argmin(np.array(distances))])

        sourceinfo['proto_dist'] = protodists
        sourceinfo['disk_dist']  = diskdists

        sourceinfo.to_csv(sourceinfofile.split('.txt')[0]+'_CoordFix.txt',sep=' ',index=False)
        
        #
        # Now the YSO Compare File
        # 

        if os.path.exists(YSOcomparefile):
            YSOinfo = pd.read_csv(YSOcomparefile,delim_whitespace=True)
            YSOinfo['RA']          = sourceinfo['RA']
            YSOinfo['DEC']         = sourceinfo['DEC']
            YSOinfo['Proto_Dist']  = protodists
            YSOinfo['Proto_Class'] = diskdists
            YSOinfo['Disk_Dist']   = protoclasses
            YSOinfo['Disk_Class']  = diskclasses

            YSOinfo.to_csv(YSOcomparefile.split('.txt')[0]+'_CoordFix.txt',sep=' ',index=False)

    #
    # Now the CalFactors File
    #

    if os.path.exists(calfactorsfile):
        calfactorsinfo = pd.read_csv(calfactorsfile,delim_whitespace=True)

        # Now this is positive, since a positive DX means a decrease in RA. Pixel value goes up, RA goes down
        calfactorsinfo['dx'] = calfactorsinfo['dx']+dx
        calfactorsinfo['dy'] = calfactorsinfo['dy']+dy

        # Don't update dx_check and dy_check 

        calfactorsinfo.to_csv(calfactorsfile.split('.txt')[0]+'_CoordFix.txt',sep=' ',index=False,float_format='%.4f')

    #
    # Now the variables file using regular expressions!
    #
    if os.path.exists(variablesfile):
        with open(variablesfile,'r') as f:
            content=f.read()
            # Find each occurence of an RA,DEC string
            each_occurence = re.findall(r'\d+.\d+,\d+.\d+',content)
            for ind,eachone in enumerate(each_occurence):
                # Again, a positive dx means the RA value should go *DOWN* (i.e. to the west)
                newRA = float(eachone.split(',')[0])-dx/3600.0
                newDEC= float(eachone.split(',')[1])+dy/3600.0
                if ind == 0:
                    content_new = re.sub(eachone,'{},{}'.format(newRA,newDEC),content)
                else:
                    content_new = re.sub(eachone,'{},{}'.format(newRA,newDEC),content_new)

        if len(each_occurence)>0: 
            with open(variablesfile.split('.txt')[0]+'_CoordFix.txt','w') as outfile:
                outfile.write(content_new)
        else:
            with open(variablesfile.split('.txt')[0]+'_CoordFix.txt','w') as outfile:
                outfile.write(content)
