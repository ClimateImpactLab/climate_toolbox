

def main(
        run_location=None,
        n_jobs=None,
        verbose=None,
        clim=None,
        shapefile_location=None,
        shapefile_name=None,
        shp_id=None,
        numeric_id_fields=None,
        string_id_fields=None,
        weightlist=None,
        use_existing_segment_shp=None,
        filter_ocean_pixels=None,
        keep_features=None,
        drop_features=None,
        run_spec_file=None):
    '''
    run_location:
        run_location
    n_jobs:
        n_jobs
    verbose:
        verbose
    clim:
        clim
    shapefile_location:
        shapefile_location
    shapefile_name:
        shapefile_name
    shp_id:
        shp_id
    numeric_id_fields:
        numeric_id_fields
    string_id_fields:
        string_id_fields
    weightlist:
        weightlist
    use_existing_segment_shp: 
        use_existing_segment_shp
    filter_ocean_pixels:
        filter_ocean_pixels
    keep_features:
        keep_features
    drop_features:
        drop_features
    run_spec_file:
        run_spec_file
    '''

    if run_spec_file is not None:
        with open(run_spec_file, 'r') as f:
            run_spec = dict(ast.literal_eval(f.read()))

    else:
        run_spec = {}


    pass
