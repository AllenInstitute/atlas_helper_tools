import urllib, json
import os

def download_metadata( dataset_id, output_path, product_id = 3 ) :

    # RMA query to get information for dataset metadata
    query_url = \
    "https://developingmouse.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,structure_sets[id$eq22],pipe::list[xstructures$eq'id'],model::SectionDataSet[id$eq%d],rma::include,genes,plane_of_section,treatments,specimen(donor(age,organism)),probes(orientation,predicted_sequence,forward_primer_sequence,reverse_primer_sequence),products[id$eq%d],model::StructureUnionize,rma::criteria,section_data_set[id$eq%d],rma::include,structure[id$in$xstructures],rma::options[only$eqid,section_data_set_id,name,expression_energy,acronym,red,green,blue],model::SectionImage[data_set_id$eq%d],rma::include,associates,alternate_images,rma::options[order$eq'sub_images.section_number$asc'],"\
    % (dataset_id,product_id,dataset_id,dataset_id)

    response = urllib.request.urlopen(query_url)
    metadata = json.loads(response.read())['msg']
    
    with open( output_path, 'w') as mfile :
        json.dump( metadata, mfile, indent=4 )

        
        
def download_images( dataset_id, output_directory, downsample_factor = 3 ) :
    
    # RMA query to get information for images in the dataset
    query_url  = "http://api.brain-map.org/api/v2/data/query.json?"
    query_url += "criteria=model::SectionImage"
    query_url += ",rma::criteria,[data_set_id$eq%d]" % (dataset_id)
    query_url += ",rma::options[num_rows$eqall]" 

    response = urllib.request.urlopen(query_url)
    images = json.loads(response.read())['msg']

    for i in images :
        #print( i['section_number'], i['id'] )
        image_url  = "http://api.brain-map.org/api/v2/section_image_download/%d?downsample=%d" % (i['id'],downsample_factor)
        image_path = os.path.join( output_directory, '%04d_%d.jpg' % (i['section_number'],dataset_id)    )
        urllib.request.urlretrieve(image_url, image_path)