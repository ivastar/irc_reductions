import threedhst
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
    
For the threedhst library:
https://code.google.com/p/threedhst/

Running Ureka ipython.

"""
    
def check_dq_irc0222b():    
    """
    Check individual FLTs and mark any blemishes with polygons.
    """
    
    import threedhst.dq
    
    files = glob.glob('IRC0222B-*asn.fits')
    for asn in files:
        threedhst.dq.checkDQ(asn_direct_file=asn, asn_grism_file=asn, path_to_flt='../RAW/')    
    
def check_dq_irc0222a():
    """
    Check individual FLTs and mark any blemishes with polygons.
    """
    
    import threedhst.dq
    
    files = glob.glob('IRC0222A-*asn.fits')
    for asn in files:
        threedhst.dq.checkDQ(asn_direct_file=asn, asn_grism_file=asn, path_to_flt='../RAW/')    

def process_images_irc0222a():
    
    import glob
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair    
    
    ### Reduce the F160W image, will use as reference, align to itself.
    file = 'IRC0222A-04-258-F160W_asn.fits'
    pair(direct_asn=file, grism_asn=None, radec=None, raw_path='../RAW/', mask_grow=8, scattered_light=False, 
    final_scale=0.06, skip_direct=False, ACS=False)
    
    for filter in ['F105W','F125W']:
        files= glob.glob('IRC0222B*'+filter+'_asn.fits')
        for i in range(len(files)):
            pair(direct_asn=files[i], grism_asn=None, radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
            skip_direct=False, ACS=False,align_threshold=8.)
    
    ### Make a catalog based on F160W image
    s = sx.SExtractor()
    s.aXeParams()
    s.sextractImage('IRC0222A-04-258-F160W_drz_sci.fits')
    sx.sexcatRegions('test.cat', 'test.reg', format=1)
    tmp_cat = sx.mySexCat('test.cat')
    radec_cat = 'IRC0222A-04-258-F160W_radec.cat'
    with open(redec_cat,'w') as f:
        for i in range(tmp_cat.nrows):
            f.write('{}\t{}\n'.format(tmp_cat['X_WORLD'][i],tmp_cat['Y_WORLD'][i]))
        
    files= glob.glob('IRC0222A-04-258-F125W_asn.fits')
    pair(direct_asn=files[i], grism_asn=None, radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
        skip_direct=False, ACS=False,align_threshold=8.)
    
    direct = glob.glob('IRC0222A*F105W_asn.fits')
    grism = glob.glob('IRC0222B*G102_asn.fits')
    for i in range(len(direct)):
        pair(direct_asn=direct[i], grism_asn=grism[i], radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
        skip_direct=False, ACS=False, align_threshold=8.)


def process_images_irc0222b():
    
    """
    Processing all images, direct and grism, for the IRC0222B cluster.
    """
    
    import glob
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair    
    import threedhst.sex as sx
        
    ### Reduce the F160W image, will use as reference, align to itself.
    file = 'IRC0222B-12-244-F160W_asn.fits'
    pair(direct_asn=file, grism_asn=None, radec=None, raw_path='../RAW/', mask_grow=8, scattered_light=False, 
    final_scale=0.06, skip_direct=False, ACS=False)
    
    ### Make a catalog based on F160W image
    s = sx.SExtractor()
    s.aXeParams()
    s.sextractImage('IRC0222B-12-244-F160W_drz_sci.fits')
    sx.sexcatRegions('test.cat', 'test.reg', format=1)
    tmp_cat = sx.mySexCat('test.cat')
    radec_cat = 'IRC0222B-12-244-F160W_radec.cat'
    with open(redec_cat,'w') as f:
        for i in range(tmp_cat.nrows):
            f.write('{}\t{}\n'.format(tmp_cat['X_WORLD'][i],tmp_cat['Y_WORLD'][i]))
    
    for filter in ['F105W','F125W']:
        files= glob.glob('IRC0222B*'+filter+'_asn.fits')
        for i in range(len(files)):
            pair(direct_asn=files[i], grism_asn=None, radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
            skip_direct=False, ACS=False,align_threshold=8.)
    
    direct = glob.glob('IRC0222B*F140W_asn.fits')
    grism = glob.glob('IRC0222B*G141_asn.fits')
    for i in range(len(direct)):
        pair(direct_asn=direct[i], grism_asn=grism[i], radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
        skip_direct=False, ACS=False, align_threshold=8.)



