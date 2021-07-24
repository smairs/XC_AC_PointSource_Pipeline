import pickle
import glob

def family_members(PSpicklefile,date_cutoff,wave):
    region = PSpicklefile.split('/')[-1].split('_')[0]
    most_recent_sourceinfo = sorted(list(glob.glob('pointsource_results/'+region+'/*_'+wave+'_sourceinfo.txt')))[-1]
    datescan = most_recent_sourceinfo.split('_')[-4]+'_'+most_recent_sourceinfo.split('_')[-3]
    outfile = open('pointsource_results/'+region+'/family_members_'+str(date_cutoff)+'_'+wave+'.txt','w')
    try:
        familyinfo = pickle.load(open(PSpicklefile,'rb'))
        family_members = familyinfo['family']
    except:
        family_members = ['NO CAL POSSIBLE']
    family_to_print = []
    for i in family_members:
        if i != 'NO CAL POSSIBLE':
            family_to_print.append(int(i))
        else:
            family_to_print.append(i)
    outfile.write('Family members derived from all data obtained before '+str(date_cutoff)+':\n\n')
    for eachmember in sorted(family_to_print):
        if eachmember != 'NO CAL POSSIBLE':
            outfile.write(str(eachmember).zfill(2)+'\n')
        else:
            outfile.write(eachmember+'\n')
    outfile.write('\n\nAll source indices correspond to the source indices presented in '+most_recent_sourceinfo)
    outfile.close()

    
