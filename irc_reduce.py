import threedhst
import glob
import os
import numpy as np
    

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
    import threedhst.sex as sx

    
    ### Reduce the F160W image, will use as reference, align to itself.
    file = 'IRC0222A-04-258-F160W_asn.fits'
    pair(direct_asn=file, grism_asn=None, radec=None, raw_path='../RAW/', mask_grow=8, scattered_light=False, 
    final_scale=0.06, skip_direct=False, ACS=False)
        
    ### Make a catalog based on F160W image
    s = sx.SExtractor()
    s.aXeParams()
    s.sextractImage('IRC0222A-04-258-F160W_drz_sci.fits')
    sx.sexcatRegions('test.cat', 'test.reg', format=1)
    tmp_cat = sx.mySexCat('test.cat')
    radec_cat = 'IRC0222A-04-258-F160W_radec.cat'
    with open(radec_cat,'w') as f:
        for i in range(tmp_cat.nrows):
            f.write('{}\t{}\n'.format(tmp_cat['X_WORLD'][i],tmp_cat['Y_WORLD'][i]))
        
    for file in glob.glob('IRC0222A-*-*-F*W_asn.fits'):    
        pair(direct_asn=file, grism_asn=None, radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
            skip_direct=False, ACS=False,align_threshold=6.)
    
    direct = ['IRC0222A-09-266-F105W_asn.fits','IRC0222A-13-256-F105W_asn.fits']
    grism = ['IRC0222A-09-266-G102_asn.fits', 'IRC0222A-13-256-G102_asn.fits']
    for i in range(len(direct)):
        pair(direct_asn=direct[i], grism_asn=grism[i], radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
        skip_direct=True, ACS=False, align_threshold=6.)
        
def make_mosaic_irc0222a():
    
    import drizzlepac
    import astropy.io.fits as fits
    import numpy as np

    ZPs = {'F105W':26.2687, 'F125W':26.25, 'F140W':26.46, 'F160W':25.96}
    
    direct_files = glob.glob('IRC0222A-*-*-F*W_asn.fits')
    f105w_files = glob.glob('IRC0222A-*-*-F105W_asn.fits')
    
    ### make one asn file with all
    threedhst.utils.combine_asn_shifts(direct_files, out_root='IRC0222A_direct',
        path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f105w_files, out_root='IRC0222A-F105W',
        path_to_FLT='./', run_multidrizzle=False)
        
    ### run astrodrizzle with all images to figure out size for mosaic
    drizzlepac.astrodrizzle.AstroDrizzle('IRC0222A_direct_asn.fits', clean=True, final_scale=0.06, 
        final_pixfrac=0.8, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', 
        driz_cr_scale = '2.5 0.7', final_wht_type = 'IVM', skysub = False, final_wcs=True)
    ### run astrodrizzle for each to make the same size mosaic
    for file in ['IRC0222A-F105W_asn.fits','IRC0222A-04-258-F125W_asn.fits','IRC0222A-04-258-F160W_asn.fits']:
        drizzlepac.astrodrizzle.AstroDrizzle(file, clean=True, context=False, final_pixfrac=0.8, preserve=False, 
            driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7', final_refimage='IRC0222A_direct_drz_sci.fits', final_wcs=True)
    
    ### coadd all images, normalizing zeropoint to F105W, weigh by inverse variance
    ### make a detection noise equalized image image
    print 'IRC0222A-F105W'
    sci_sum = fits.open('IRC0222A-F105W_drz_sci.fits')
    wht_sum = fits.open('IRC0222A-F105W_drz_wht.fits')
    
    sci_sum[0].data = sci_sum[0].data*wht_sum[0].data
    noise_eq = sci_sum[0].data*np.sqrt(wht_sum[0].data)
    
    for root, filter in zip(['IRC0222A-04-258-F125W','IRC0222A-04-258-F160W'],['F125W','F160W']):
        print root
        sci = fits.open('{}_drz_sci.fits'.format(root))
        wht = fits.open('{}_drz_wht.fits'.format(root))
        index = sci[0].data > 1.e6
        sci[0].data[index] = 0
        wht[0].data[index] = 0
        
        zp_factor = 10**((ZPs['F105W']-ZPs[filter])/2.5)
        sci[0].data = sci[0].data*zp_factor
        wht[0].data = wht[0].data/zp_factor**2
        
        sci_sum[0].data += sci[0].data*wht[0].data
        wht_sum[0].data += wht[0].data
        noise_eq += sci[0].data*np.sqrt(wht[0].data)
        
        sci.close()
        wht.close()
        
    del(sci)
    del(wht)
        
    index = wht_sum[0].data == 0
    sci_full = sci_sum[0].data/wht_sum[0].data
    sci_full[index] = 0

    print 'Writing final images.'
    fits.writeto('IRC0222A-IR_sci.fits', data=sci_full, header=sci_sum[0].header, clobber=True)
    fits.writeto('IRC0222A-IR_wht.fits', data=wht_sum[0].data, header=sci_sum[0].header, clobber=True)
    fits.writeto('IRC0222A_noise_equalized.fits', data=noise_eq, header=sci_sum[0].header, clobber=True)
        
    ### run sextractor to make caalog

    sextr = "sex %s -c %s.config -WEIGHT_IMAGE %s" %('IRC0222A-IR_sci.fits','GOODS-S_F160W_v1','IRC0222A-IR_wht.fits')
    os.system(sextr)
    

def copy_flts(field='IRC0222A'):
    
    if field == 'IRC0222A':
        files = ['IRC0222A-09-266-F105W_asn.fits','IRC0222A-13-256-F105W_asn.fits','IRC0222A-09-266-G102_asn.fits','IRC0222A-13-256-G102_asn.fits']
    if field == 'IRC0222B':
        files = ['IRC0222B-05-244-F140W_asn.fits','IRC0222B-10-254-F140W_asn.fits','IRC0222B-05-244-G141_asn.fits','IRC0222B-10-254-G141_asn.fits']
    
    for asn_file in files:
        os.system('rsync -av {} ../INTERLACE_{}/'.format(asn_file,field))
        asn = threedhst.utils.ASNFile(asn_file)
        for exp in asn.exposures:
            os.system('rsync -av {}_flt.fits ../INTERLACE_{}/'.format(exp,field))


def interlace_irc0222a():
    """
    Interlace the final FLT images and make an interlaced reference image.
    Create a model. Refine background.
    Extract objects down to F105W=24.
    This is all done in the INTERLACE_IRC0222A directory.
    """

    import unicorn
    from unicorn.reduce import adriz_blot_from_reference as adriz_blot
    import scipy.ndimage as nd
    from astropy.io import fits
             
    NGROWX=100
    NGROWY=1
    pad=60
    CATALOG='../PREP_FLT/sextr/IRC0222A-IR.cat'
    REF_IMAGE = '../PREP_FLT/IRC0222A-IR_sci.fits'
    SEG_IMAGE = '../PREP_FLT/sextr/IRC0222A-IR.seg.fits'
    REF_FILTER='F105W'
    REF_EXT = 0

    grism=glob.glob('IRC0222A-*-*-G102_asn.fits')
            
    extract_limit = 35.
    skip_completed=False

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(grism)):
        pointing=grism[i].split('-G102')[0]
        adriz_blot(pointing=pointing+'-F105W', pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, growx=2, growy=2, auto_offsets=True, ref_exp=0, ref_image=REF_IMAGE, ref_ext=REF_EXT, ref_filter=REF_FILTER, seg_image=SEG_IMAGE, cat_file=CATALOG, grism='G102')
        unicorn.reduce.interlace_combine(pointing+'-F105W', view=False, use_error=True, make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True, ref_exp=0)
        unicorn.reduce.interlace_combine(pointing+'-G102', view=False, use_error=True, make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True, ref_exp=0)
    
    # Make models.
    inter = glob.glob('IRC0222A-*G102_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G102_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=26., REFINE_MAG_LIMIT = 23.,
                make_zeroth_model=False, BEAMS=['A','B','C','D','E'],
                grow_factor=2, growx=2, growy=2, direct='F105W', grism='G102')
            if not os.path.exists(os.path.basename(model.root) + '-G102_maskbg.dat'):
                 model.refine_mask_background(grow_mask=12, threshold=0.001, update=True,
                     resid_threshold=4, clip_left=640, save_figure=True, interlace=True)
       
    # Extract objects.
    inter = glob.glob('IRC0222A-*G102_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G102_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, grism='G102',direct='F105W')
        model.extract_spectra_and_diagnostics(MAG_LIMIT=24.)
    
    
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
    with open(radec_cat,'w') as f:
        for i in range(tmp_cat.nrows):
            f.write('{}\t{}\n'.format(tmp_cat['X_WORLD'][i],tmp_cat['Y_WORLD'][i]))
    
    for filter in ['F105W','F125W','F140W','F160W']:
        files= glob.glob('IRC0222B*'+filter+'_asn.fits')
        print files
        for i in range(len(files)):
            pair(direct_asn=files[i], grism_asn=None, radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
            skip_direct=False, ACS=False,align_threshold=6.)
    
    direct = glob.glob('IRC0222B*F140W_asn.fits')
    grism = glob.glob('IRC0222B*G141_asn.fits')
    for i in range(len(direct)):
        pair(direct_asn=direct[i], grism_asn=grism[i], radec=radec_cat, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, 
        skip_direct=False, ACS=False, align_threshold=8.)


def make_mosaic_irc0222b():
    
    import drizzlepac
    import astropy.io.fits as fits
    import numpy as np

    ZPs = {'F105W':26.2687, 'F125W':26.25, 'F140W':26.46, 'F160W':25.96}
    
    direct_files = glob.glob('IRC0222B-*-*-F*W_asn.fits')
    f140w_files = glob.glob('IRC0222B-*-*-F140W_asn.fits')
    
    ### make one asn file with all
    threedhst.utils.combine_asn_shifts(direct_files, out_root='IRC0222B_direct',
        path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f140w_files, out_root='IRC0222B-F140W',
        path_to_FLT='./', run_multidrizzle=False)
    
    #def combine_images():
    ### run astrodrizzle with all images to figure out size for mosaic
    drizzlepac.astrodrizzle.AstroDrizzle('IRC0222B_direct_asn.fits', clean=True, final_scale=0.06, 
        final_pixfrac=0.8, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', 
        driz_cr_scale = '2.5 0.7', final_wht_type = 'IVM', skysub = False, final_wcs=True)
    ### run astrodrizzle for each to make the same size mosaic
    for file in ['IRC0222B-11-244-F105W_asn.fits','IRC0222B-12-244-F125W_asn.fits','IRC0222B-F140W_asn.fits','IRC0222B-12-244-F160W_asn.fits']:
        drizzlepac.astrodrizzle.AstroDrizzle(file, clean=True, context=False, final_pixfrac=0.8, preserve=False, 
            driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7', final_refimage='IRC0222B_direct_drz_sci.fits', final_wcs=True)
    
    
    ### coadd all images, normalizing zeropoint to F105W, weigh by inverse variance
    ### make a detection noise equalized image image
    print 'IRC0222B-F140W'
    sci_sum = fits.open('IRC0222B-F140W_drz_sci.fits')
    wht_sum = fits.open('IRC0222B-F140W_drz_wht.fits')
    
    sci_sum[0].data = sci_sum[0].data*wht_sum[0].data
    noise_eq = sci_sum[0].data*np.sqrt(wht_sum[0].data)
    
    for root, filter in zip(['IRC0222B-11-244-F105W','IRC0222B-12-244-F125W','IRC0222B-12-244-F160W'],['F105W','F125W','F160W']):
        print root
        sci = fits.open('{}_drz_sci.fits'.format(root))
        wht = fits.open('{}_drz_wht.fits'.format(root))
        index = sci[0].data > 1.e6
        sci[0].data[index] = 0
        wht[0].data[index] = 0
        
        zp_factor = 10**((ZPs['F105W']-ZPs[filter])/2.5)
        sci[0].data = sci[0].data*zp_factor
        wht[0].data = wht[0].data/zp_factor**2
        
        sci_sum[0].data += sci[0].data*wht[0].data
        wht_sum[0].data += wht[0].data
        noise_eq += sci[0].data*np.sqrt(wht[0].data)
        
        sci.close()
        wht.close()
        
    del(sci)
    del(wht)
        
    index = wht_sum[0].data == 0
    sci_full = sci_sum[0].data/wht_sum[0].data
    sci_full[index] = 0

    print 'Writing final images.'
    fits.writeto('IRC0222B-IR_sci.fits', data=sci_full, header=sci_sum[0].header, clobber=True)
    fits.writeto('IRC0222B-IR_wht.fits', data=wht_sum[0].data, header=sci_sum[0].header, clobber=True)
    fits.writeto('IRC0222B_noise_equalized.fits', data=noise_eq, header=sci_sum[0].header, clobber=True)
        
    ### run sextractor to make caalog

    sextr = "sex %s -c %s.config -WEIGHT_IMAGE %s" %('IRC0222B-IR_sci.fits','IRC0222B','IRC0222B-IR_wht.fits')
    os.system(sextr)

def interlace_irc0222b():

    """
    Interlace the final FLT images and make an interlaced reference image.
    Create a model. Refine background.
    Extract objects down to F140W=24.
    This is all done in the INTERLACE_IRC0222B directory.
    """
    
    import unicorn
    from unicorn.reduce import adriz_blot_from_reference as adriz_blot
    import scipy.ndimage as nd
    from astropy.io import fits
             
    NGROWX=100
    NGROWY=1
    pad=60
    CATALOG='../PREP_FLT/sextr/IRC0222B-IR.cat'
    REF_IMAGE = '../PREP_FLT/IRC0222B-IR_sci.fits'
    SEG_IMAGE = '../PREP_FLT/sextr/IRC0222B-IR.seg.fits'
    REF_FILTER='F140W'
    REF_EXT = 0

    grism=glob.glob('IRC0222B-*-*-G141_asn.fits')
            
    extract_limit = 35.
    skip_completed=False

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(grism)):
        pointing=grism[i].split('-G141')[0]
        adriz_blot(pointing=pointing+'-F140W', pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, growx=2, growy=2, auto_offsets=True, ref_exp=0, ref_image=REF_IMAGE, ref_ext=REF_EXT, ref_filter=REF_FILTER, seg_image=SEG_IMAGE, cat_file=CATALOG, grism='G141')
        unicorn.reduce.interlace_combine(pointing+'-F140W', view=False, use_error=True, make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True, ref_exp=0)
        unicorn.reduce.interlace_combine(pointing+'-G141', view=False, use_error=True, make_undistorted=False, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True, ref_exp=0)
    
    #### Make models.
    inter = glob.glob('IRC0222B-*G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=26., REFINE_MAG_LIMIT = 23.,
                make_zeroth_model=False, BEAMS=['A','B','C','D','E'],
                grow_factor=2, growx=2, growy=2, direct='F140W', grism='G141')
            if not os.path.exists(os.path.basename(model.root) + '-G141_maskbg.dat'):
                 model.refine_mask_background(grow_mask=12, threshold=0.001, update=True,
                     resid_threshold=4, clip_left=640, save_figure=True, interlace=True)

    #### Extract objects.
    inter = glob.glob('IRC0222B-*G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, grism='G141')
            model.extract_spectra_and_diagnostics(MAG_LIMIT=24.)
