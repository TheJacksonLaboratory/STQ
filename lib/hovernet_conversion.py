import os
import json
import pandas as pd
import numpy as np

def getSpotCellData(id, downSamplingFactor=1.13215, factor=1, space_ranger_output_path='', hovernet_output_path='', data_save_path='',
                   classes={0: 'no_label', 1: 'neoplasia', 2: 'inflammatory', 3: 'connective', 4: 'necrosis', 5: 'no_neoplasia'}):
    
    """downSamplingFactor: some factor that was given (?) to hovernet to reduce image dimensions. Need it here to rescale the coordinates
    and match Space Ranger alignment output files
    
    factor: if 1 then cells under ST spots are considered, multiplies the spot diameter
    
    space_ranger_output_path: directoriers with two files for each id: /<id>/spatial/scalefactors_json.json and /spatial/tissue_positions_list.csv
    
    hovernet_output_path: directory with files with extract_hov_counts.csv for each id
    
    data_save_path: directory to save output to
    """
    
    print('Processing ID: %s' % id)
    
    df_spaceranger_alignment = pd.read_csv(space_ranger_output_path + '%s/spatial/tissue_positions_list.csv' % id, index_col=0, header=None)
    # Format description at: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    df_spaceranger_alignment.columns = ['in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
    df_spaceranger_alignment.index.name = 'barcode'

    scalefactors_json_file = space_ranger_output_path + '%s/spatial/scalefactors_json.json' % id
    with open(scalefactors_json_file) as tempfile:
        scalefactors_tbl = json.load(tempfile)

    files = os.listdir(hovernet_output_path)
    hfile = [file for file in files if ((id in file) and ('.extract_hov_counts.csv' in file))][0]

    hovernet_cells = pd.read_csv(hovernet_output_path + hfile, index_col=0)
    hovernet_cells['cell_type'] = hovernet_cells['cell_type'].replace(classes)

    def get_st_barcode(centroid_y, centroid_x, df_spaceranger_alignment, spot_diameter_fullres):

        centroid_x *= downSamplingFactor
        centroid_y *= downSamplingFactor

        se_i = (df_spaceranger_alignment['pxl_row_in_fullres'] - centroid_x).abs() < (spot_diameter_fullres/2)
        se_j = (df_spaceranger_alignment['pxl_col_in_fullres'] - centroid_y).abs() < (spot_diameter_fullres/2)
        se = se_i & se_j
        matched = se[se].index

        if len(matched)>1:
            print('Decrease ovelap factor to have unique barcode matching')
            print(matched)
        assert len(matched)<=1

        if len(matched)>0:
            barcode = matched[0]
        else:
            barcode = None

        return barcode

    suffFactor = '' if factor==1 else 'factor_%s_' % factor

    hovernet_cells['barcode'] = hovernet_cells[['centroid_x', 'centroid_y']].apply(lambda coords: get_st_barcode(*coords, df_spaceranger_alignment, factor * scalefactors_tbl['spot_diameter_fullres']), axis=1)
    hovernet_cells.to_csv(data_save_path + hfile + '.%swith_barcode.csv.gz' % suffFactor)
    
    df_filtered = hovernet_cells.dropna().reset_index().set_index('barcode')

    df_class_densities = pd.get_dummies(df_filtered['cell_type']).groupby(level=0).mean()
    df_class_densities.columns = 'density_' + pd.Series(df_class_densities.columns).replace(classes)

    df_class_counts = pd.get_dummies(df_filtered['cell_type']).replace({0: np.nan}).groupby(level=0).count()
    df_class_counts.columns = 'count_' + pd.Series(df_class_counts.columns).replace(classes)

    df_total_counts = df_class_counts.sum(axis=1).to_frame()
    df_total_counts.columns = ['total_count']

    df_area_classes = []
    for barcode in df_filtered.index.unique():
        df_temp = df_filtered.loc[barcode]
        if type(df_temp) is pd.Series:
            df_temp = df_temp.to_frame().T
        df_area_classes.append(df_temp.set_index('cell_type')['cell_type_area'].groupby(level=0).mean().rename(barcode))
    df_area_classes = pd.concat(df_area_classes, axis=1).T
    df_area_classes.columns = 'average_area_' + df_area_classes.columns
    df_area_classes.index.name = 'barcode'

    df_area = df_filtered['cell_type_area'].groupby(level=0).mean().to_frame()
    df_area.columns = ['average_area']

    df_added = pd.concat([df_class_densities, df_class_counts, df_total_counts, df_area_classes, df_area], axis=1)

    df_added.to_csv(data_save_path + '%s_%shovenet_info_per_barcode.csv' % (id, suffFactor))

    return
