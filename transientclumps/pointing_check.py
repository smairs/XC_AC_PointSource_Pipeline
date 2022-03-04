from transientclumps.TCGaussclumpsFunctions import *
from transientclumps.TCOffsetFunctions import *
from transientclumps.TCPrepFunctions import *
from transientclumps.merge_catalog import *
import glob as glob
import matplotlib.pyplot as plt

def pointing_check(direc,wave):

    sdf_files = sorted(glob.glob(direc+'/*'+wave+'*mJybmsm.sdf'))
    print('\n\n\n',sdf_files,'\n\n\n')

    if wave == '450':
        pix_scale = 2.0
    elif wave == '850':
        pix_scale = 3.0

    region = sdf_files[0].split('/')[-1].split('_')[0]

    xoffs = []
    yoffs = []
    peakrs = []
    xoffuncs = []
    yoffuncs = []
    peakruncs = []
    datescans = []

    for eachfile in sdf_files:

        #eachfile = eachfile

        datescans.append(eachfile.split('_00')[0].split('_')[-1]+'_00'+eachfile.split('_00')[-1].split('_')[0])

        prepare_image(eachfile)

        eachfile = eachfile.split('/')[-1]

        run_gaussclumps(eachfile.split('_mJybmsm')[0]+'_Jypbmsm_crop.sdf','parameters/GCParms.txt')

        target_catalogue = eachfile.split('_mJybmsm')[0]+'_Jypbmsm_crop'+'_log.FIT'
        reference_catalogue = sdf_files[0].split('/')[-1].split('_mJybmsm')[0]+'_Jypbmsm_crop'+'_log.FIT'

        offsets, errs, NSurviveCull, Nmatch, MatchInd, TargIndMatched, CulledInd, TargIndCulled  = source_match(target_catalogue, reference_catalogue, minpeak=0.2, maxrad=30, maxsep=10, cutoff=4, pix_scale=pix_scale)

        avg_xoff  = offsets[0]
        avg_yoff  = offsets[1]
        avg_peakr = offsets[2]
        std_xoff  = errs[0]
        std_yoff  = errs[1]
        std_peakr = errs[2]

        xoffs.append(avg_xoff)
        yoffs.append(avg_yoff)
        peakrs.append(avg_peakr)
        xoffuncs.append(std_xoff)
        yoffuncs.append(std_yoff)
        peakruncs.append(std_peakr)


        cat_match_name = target_catalogue.split('.FIT')[0]+'_match.FIT'

        #print(eachfile,target_catalogue,reference_catalogue,eachfile.split('_00')[0].split('_')[-1],eachfile.split('_00')[-1].split('_')[0])

        merge_catalog(target_catalogue, reference_catalogue, eachfile.split('_00')[0].split('_')[-1],str(int(eachfile.split('_00')[-1].split('_')[0])), cat_match_name,ref_index=MatchInd,cat_index=TargIndMatched)

        cat_cull_name = target_catalogue.split('.FIT')[0]+'_cull.FIT'
        
        merge_catalog(target_catalogue, reference_catalogue, eachfile.split('_00')[0].split('_')[-1],str(int(eachfile.split('_00')[-1].split('_')[0])), cat_cull_name,ref_index=CulledInd,cat_index=TargIndCulled)

    prev = False
    if os.path.exists('pointsource_results/'+region+'/'+region+'_'+wave+'_local_pointing_corrections.txt'):
        prev = True
        previous_file = open('pointsource_results/'+region+'/'+region+'_'+wave+'_local_pointing_corrections.txt','r')
        content = previous_file.readlines()

    pcheckfile = open('pointsource_results/'+region+'/'+region+'_'+wave+'_local_pointing_corrections.txt','w')

    all_xoffs = []
    all_yoffs = []
    all_xoffuncs = []
    all_yoffuncs = []

    if prev:
        for i in content:
            pcheckfile.write(i)
            if i[0] != '#':
                all_xoffs.append(float(i.split('\t')[1]))
                all_yoffs.append(float(i.split('\t')[2]))
                all_xoffuncs.append(float(i.split('\t')[3]))
                all_yoffuncs.append(float(i.split('\t')[4].split('\n')[0]))

    else:
        pcheckfile.write('#Datescan\tXoff\tYoff\tXoffunc\tYoffunc\n')

    dummy = -1
    for i,j,k,l,m in zip(datescans,xoffs,yoffs,xoffuncs,yoffuncs):
        dummy+=1
        print(i,j,k,l,m)
        if dummy == 0:
            if prev:
                continue # This ensures we don't add the first (0,0) point multiple times to the lines to the pointing check file
        pcheckfile.write(str(i)+'\t'+str(j)+'\t'+str(k)+'\t'+str(l)+'\t'+str(m)+'\n')
        all_xoffs.append(float(j))
        all_yoffs.append(float(k))
        all_xoffuncs.append(float(l))
        all_yoffuncs.append(float(m))

    pcheckfile.close()
        

    os.system('mv *_Jypbmsm_crop* '+direc)
    plt.errorbar(all_xoffs,all_yoffs,xerr=all_xoffuncs,yerr=all_yoffuncs,linestyle='None',marker='o')
    plt.xlabel('RA Offsets (arcsecs)')
    plt.ylabel('Dec Offsets (arcsecs)')
    plt.savefig('pointsource_results/'+region+'/'+region+'_'+wave+'_local_pointing_corrections.pdf',format='pdf')

    plt.clf()

    #plt.hist(peakrs)
    #plt.xlabel('Flux Ratio with Reference')
    #plt.ylabel('Count')
    #plt.savefig('pointsource_results/'+region+'/'+region+'_'+wave+'_local_flux_ratio_with_ref.pdf',format='pdf')

    #plt.clf()
