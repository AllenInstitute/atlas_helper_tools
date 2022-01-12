import json
import numpy as np
import SimpleITK as sitk
import os

REDCHANNEL = 0
GREENCHANNEL = 1
BLUECHANNEL = 2

#
# References:
#
# http://help.brain-map.org/display/mousebrain/API
# http://help.brain-map.org/download/attachments/2818169/ABADataProductionProcesses.pdf
# http://help.brain-map.org/download/attachments/2818169/InformaticsDataProcessing.pdf
#
# http://help.brain-map.org/display/devmouse/API
# https://help.brain-map.org/download/attachments/4325389/DevMouse_Overview.pdf
# https://help.brain-map.org/download/attachments/4325389/DevMouse_InformaticsDataProcessing.pdf
#

def read_atlas( atlas_directory, age, plane_of_section ) :
    
    # References:
    # http://help.brain-map.org/display/devmouse/API#API-DownloadAtlas3-DReferenceModels
    # http://help.brain-map.org/display/mousebrain/API#API-DownloadAtlas3-DReferenceModels
    
    image_file = os.path.join( atlas_directory, age, plane_of_section, 'atlasVolume.nii.gz' )
    volume = sitk.ReadImage( image_file )
    
    xform_file = os.path.join( atlas_directory, age, plane_of_section, 'CanonicalToAtlasVolume.tfm' )
    xform = sitk.ReadTransform( xform_file )
    
    atlas = {}
    atlas['volume'] = volume
    atlas['canonical_to_atlas_transform'] = xform
    
    return atlas

def read_alignment3d( json ) :
    
    # References:
    # http://api.brain-map.org/doc/Alignment3d.html
    # https://itk.org/Doxygen/html/classitk_1_1AffineTransform.html
    
    parameters = {}
    for t in ['trv','tvr'] :
        keys = ['%s_%02d' % (t,p) for p in range(12)]
        parameters[t] = [json[k] for k in keys]    
    return parameters

def read_alignment2d( json ) :
    
    # References:
    # http://api.brain-map.org/doc/Alignment2d.html
    # https://itk.org/Doxygen/html/classitk_1_1AffineTransform.html
    
    parameters = {}
    for t in ['tsv','tvs'] :
        keys = ['%s_%02d' % (t,p) for p in range(6)]
        parameters[t] = [json[k] for k in keys]    
    return parameters
    
def compute_section_info( images ) :
    
    indices = [x['section_number'] for x in images]
    indices = np.sort(indices)
    min_section_number = indices[0]
    section_span = indices[-1] - indices[0]
    
    diff = np.diff(indices)
    diff = np.sort(diff)
    min_section_gap = diff[0]
    
    resolutions = [x['resolution'] for x in images]
    min_resolution = np.unique( resolutions )
    if len(min_resolution) :
        "Warning: images have different resolutions - using smallest value"
    min_resolution = min_resolution[0]
    
    return( section_span, min_section_number, min_section_gap, min_resolution )
    
    
def read_json( file_path ) :
    
    with open(file_path,"r") as file:
        payload = json.load(file)
        
    payload = payload[0]
    
    data = {}
    data['id'] = payload['id']
    data['treatment'] = payload['treatments'][0]['name']
    data['plane_of_section'] = payload['plane_of_section']['name']
    data['age'] = payload['specimen']['donor']['age']['name']
    data['section_thickness'] = payload['section_thickness']
    
    data['alignment3d'] = read_alignment3d(payload['alignment3d'])
    
    images = []
    
    for s in payload['section_images'] :
    
        rec = {}
        attrs = ['section_number','width','height','resolution']
        for a in attrs :
            rec[a] = s[a]
            
        rec['alignment2d'] = read_alignment2d(s['alignment2d'])
        images.append(rec)
        
    data['section_images'] = images
    
    (ss, msn, msg, mr) = compute_section_info( images )
    data['section_span'] = ss
    data['min_section_number'] = msn
    data['min_section_gap'] = msg
    data['min_resolution'] = mr
    data['section_spacing'] = msg * data['section_thickness']
    data['section_origin'] = msn * data['section_thickness']
    
    return data


def initialize_output( data, atlas_directory, downsample_factor = 3 ) :
    
    df = pow(2,downsample_factor)
    xyres = data['min_resolution'] * df
    zres = data['section_spacing']
    
    json_file = os.path.join( atlas_directory, 'atlas_metadata.json' )
       
    with open(json_file,"r") as file:
        payload = json.load(file)
        
    ssize = payload[data['age']][data['plane_of_section']]['standard_size']
    xsize = int(round(ssize['width'] / xyres, 0 ))
    ysize = int(round(ssize['height'] / xyres, 0 ))
    zsize = int(round((data['section_span'] / data['min_section_gap']) + 1, 0))
    
    zorig = data['section_origin']
    
    image3d = sitk.Image( (xsize,ysize,zsize), sitk.sitkUInt8)
    image3d.SetSpacing( (xyres,xyres,zres) )
    image3d.SetOrigin( (0.0,0.0,zorig) )
    
    mask3d = sitk.Image( (xsize,ysize,zsize), sitk.sitkUInt8)
    mask3d.SetSpacing( (xyres,xyres,zres) )
    mask3d.SetOrigin( (0.0,0.0,zorig) )
    
    image2d = sitk.Image( image3d.GetSize()[0:2], sitk.sitkUInt8 )
    image2d.SetSpacing( image3d.GetSpacing()[0:2] )
    image2d.SetOrigin( image3d.GetSpacing()[0:2] )
        
    output = {}
    output['volume'] = image3d
    output['mask'] = mask3d
    output['reference_image'] = image2d
    
    return output

  
def populate_output( data, image_directory, downsample_factor, output, channel = BLUECHANNEL ) :
    
    blueChannelIndex = 2
    xyres = pow(2,downsample_factor)
    
    for x in data['section_images'] :
        
        basename = '%04d_%d.jpg' % (x['section_number'],data['id'])
        image_file = os.path.join( image_directory, basename )
        img = sitk.ReadImage( image_file )
        img = sitk.VectorIndexSelectionCast(img,channel,sitk.sitkUInt8)
        img = sitk.InvertIntensity(img)
        img.SetSpacing((xyres,xyres))
        
        # tvs = transform from volume slice to image pixel
        tvs = sitk.AffineTransform(2)
        tvs.SetParameters( x['alignment2d']['tvs'])
        
        resampled = sitk.Resample( img, output['reference_image'], tvs )
        
        zpoint = x['section_number'] * data['section_thickness']
        zindex = output['volume'].TransformPhysicalPointToIndex((0,0,zpoint))[2]
        
        output['volume'][:,:,zindex] = resampled
        output['mask'][:,:,zindex] = 1
        
        
def resample_to_atlas( data, atlas, output ) :
    
    # tvr = transform from volume to canonical reference
    # trv = transform from canonical reference to volume
    # tar = transform from atlas to canonical reference
    # tra = transform from canonical reference to atlas

    trv = sitk.AffineTransform(3)
    trv.SetParameters( data['alignment3d']['trv'] )
    
    tvr = sitk.AffineTransform(3)
    tvr.SetParameters( data['alignment3d']['tvr'] )
    
    tra = atlas['canonical_to_atlas_transform']
    tar = tra.GetInverse()
    
    # tva = tvr o tra = transfrom from volume to atlas
    # tav = tar o trv = transform from atlas to volume
    
    tva = sitk.CompositeTransform( tra )
    tva.AddTransform( tvr )
    
    tav = sitk.CompositeTransform( trv )
    tav.AddTransform( tar )
    
    resampled = sitk.Resample( output['volume'], atlas['volume'], tav )
    output['resampled_volume'] = resampled
    
    resampled = sitk.Resample( output['mask'], atlas['volume'], tav )
    output['resampled_mask'] = resampled
    
    affine3d = sitk.AffineTransform(3)
    affine3d.SetParameters( data['alignment3d']['tvr'] )

    resampled = sitk.Resample( atlas['volume'], output['volume'], tva )
    output['resampled_atlas'] = resampled


def write_volumes( output, output_directory ) :

    olist = ['volume','mask','resampled_volume','resampled_mask','resampled_atlas']
    for x in olist :
        output_file = os.path.join( output_directory, '%s.nii.gz' % x )
        sitk.WriteImage( output[x], output_file, True )
        
