    import threedhst
    import unicorn
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
    
    """
    unicorn.candels.make_asn_files(uniquename=True)
    The following are the ASN files:
    
    IRC0222A-09-266-F105W 5
    IRC0222A-13-256-F105W 5
    IRC0222A-04-258-F125W 8
    IRC0222A-04-258-F160W 8
    IRC0222A-09-266-G102 12
    IRC0222A-13-256-G102 12 
    IRC0222B-11-244-F105W 8
    IRC0222B-12-244-F125W 8
    IRC0222B-10-254-F140W 4
    IRC0222B-05-244-F140W 4
    IRC0222B-12-244-F160W 8
    IRC0222B-05-244-G141 8
    IRC0222B-10-254-G141 8
    
    """
    

def process_images_irc0222b():
    
    import unicorn.candels
    import glob
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair    
        
    for filter in ['F105W','F140W','F125W','F160W']:
        files= glob.glob('IRC0222B*'+filter+'_asn.fits')
        for i in range(len(files)):
            pair(direct_asn=files[i], grism_asn=None, radec=None, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
            skip_direct=False, ACS=False)
    
