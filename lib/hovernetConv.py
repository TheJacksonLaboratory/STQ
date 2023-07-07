import os
import json
import pandas as pd
import numpy as np
import cv2

from skimage.segmentation import expand_labels

def close_contour(c):

    def lfunc(s):
        l = np.abs(np.diff(s)).max()
        s = np.vstack([np.linspace(s[0, 0], s[0, 1], l+1, dtype=s.dtype),
                       np.linspace(s[1, 0], s[1, 1], l+1, dtype=s.dtype)])
        return s
        
    c = np.array(c)  
     
    return np.concatenate([lfunc(c[i:i+2].T).T[1:] for i in range(len(c)-1)] + [lfunc(np.vstack([c[-1], c[0]]).T).T])

def minimize_contour(c):

    cc = np.vstack([c, c[:2]])
    
    toKeep = []
    for i in range(len(c)):
        x_prev, y_prev = cc[i]
        x_curr, y_curr = cc[i+1]
        x_next, y_next = cc[i+2]
        
        keep = True
        if np.abs(x_prev-x_next)==np.abs(y_prev-y_next) and np.abs(x_prev-x_next)>=2:
            keep=False
        elif np.abs(x_prev-x_next)>=2 and y_prev==y_next and y_prev==y_curr:
            keep = False
        elif x_prev==x_next and x_prev==x_curr and np.abs(y_prev-y_next)>=2:
            keep = False   
                    
        toKeep.append(keep)
        
    return c[np.roll(toKeep, 1)]

def makeJSONoutput(stardist_details):

    """ {'mag': 40, 'nuc': {'1': {'centroid': [x_float, y_float], 'contour': [[x_int, y_int], ...], 'type_prob': 0.995, 'type': 0}, '2': ...}}
    """

    data = {'mag': 40, 'nuc': dict()}
    for i in range(len(stardist_details['prob'])):
        data['nuc'].update({str(i): dict()})
        data['nuc'][str(i)].update({'centroid': stardist_details['points'][i].astype(int)[::-1].tolist()})
        data['nuc'][str(i)].update({'contour':  minimize_contour(close_contour(stardist_details['coord'][i].T.astype(int)))[:,::-1].tolist()})
        data['nuc'][str(i)].update({'type_prob': np.float64(stardist_details['prob'][i])})
        data['nuc'][str(i)].update({'type': 0})
        
    return data

def checkMask(thumbPath, maskPath, savePath, bc = 210.):

    import matplotlib.pyplot as plt
    from skimage.transform import resize

    thumb = plt.imread(thumbPath)
    original_mask = plt.imread(maskPath)

    mask = resize(original_mask, (thumb.shape[0], thumb.shape[1]), order=3)
    mask[mask>=0.5] = 1.
    mask[mask<0.5] = 0.

    b = thumb.mean(axis=2)[mask==0].mean()
    print('Average background:', b)

    if b < bc:       
        original_mask[:] = 1.
    
    if not os.path.exists(savePath):
        os.makedirs(savePath)
    
    cv2.imwrite(savePath + os.path.basename(maskPath), original_mask * 255)

    return

def computeContourQuantities(c, na_filler: float = np.nan, scale_factor: float = 1.0):
    
    """Compute area and perimeter length, derive orientation and eccentricity for a contour
    
    Parameters:
    c: numpy.rray with shape (N, 2), closed nucleus contour
    
    na_filler: filler for undefined quantities, e.g. orientation is undefined when contour is a circle
    
    Output:
    Dictionary
    """
    
    results = dict()
    
    try:
        _length = cv2.arcLength(c, closed=True)
    except:
        _length = na_filler
        
    results['perimeter_length'] = _length * scale_factor
    
    # https://en.wikipedia.org/wiki/Image_moment
    # https://docs.opencv.org/3.4/d8/d23/classcv_1_1Moments.html
    moments = cv2.moments(c)
    
    results['area'] = moments['m00'] * (scale_factor**2)
    
    try:
        mu20p = moments['mu20'] / moments['m00']
        mu02p = moments['mu02'] / moments['m00']
        mu11p = moments['mu11'] / moments['m00']
    except:
        mu02p = mu20p = mu11p = na_filler
    
    try:
        orientation = 0.5 * np.arctan(2 * mu11p / (mu20p - mu02p))
        orientation *= 180. / np.pi
    except:
        orientation = na_filler
        
    results['orientation'] = orientation
    
    try:
        lp1 = (mu20p + mu02p) / 2.
        lp2 = 0.5 * np.sqrt(4 * mu11p**2 + (mu20p - mu02p)**2)
        lambda1 = lp1 + lp2
        lambda2 = lp1 - lp2
        eccentricity = np.sqrt(1 - (lambda2 / lambda1))
    except:
        eccentricity = na_filler
        
    results['eccentricity'] = eccentricity
    
    results['circularity'] =  4 * np.pi * results['area'] / (results['perimeter_length'] ** 2)
    
    return results


def loadNuclei(json_file_path: str = '', savepath: str = None, sname: str = '', original_mpp: float = 0.25):
    
    """Load JSON file produced by HoVer-Net inference.
    
    Parameters:
    json_file_path: path to JSON file that contains nuclei segmentation and classification results
        
    savepath: path to save data
        
    sname: identifier for a file to save
        
    Output:
    pandas.DataFrame
    """
    
    with open(json_file_path, 'r') as temp_file:
        json_data = json.loads(temp_file.read())
    
    #hovernet_magnification = float(json_data['mag'])
    #scale_factor = original_mpp / hovernet_mpp
    scale_factor = 1.
    
    df = pd.DataFrame({id: [json_data['nuc'][id]['centroid'],
                            json_data['nuc'][id]['type'],
                            json_data['nuc'][id]['type_prob']] for id in json_data['nuc'].keys()}).T.rename_axis('ids')
    df[['centroid_x', 'centroid_y']] = (df[0].apply(pd.Series) * scale_factor).astype(int)
    df['cell_type'] = df[1]
    df['cell_type_prob'] = df[2]
    df = df.drop([0, 1, 2], axis=1)
   
    contours = {id: json_data['nuc'][id]['contour'] for id in df.index}
    
    df_info = pd.DataFrame({id: computeContourQuantities(np.array(contours[id]), scale_factor=scale_factor) for id in list(contours)}).T.rename_axis('ids')
    
    df = pd.concat([df, df_info], axis=1)
    
    if not savepath is None:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
            
        df.to_csv(savepath + '%s.csv.gz' % sname)
    
    return df


def assignNuceiToSpots(grid_file_path: str = None, scalefactors_json_file: str = None, downSamplingFactor=1.13215, factor=1, spot_shape='circle',
                    space_ranger_output_path='', hovernet_data_file_path='', savepath='',
                   classes={0: 'no_label', 1: 'neoplasia', 2: 'inflammatory', 3: 'connective', 4: 'necrosis', 5: 'no_neoplasia'}):
    
    """
    grid_file_path: path to file 'tissue_positions_list.csv' from spaceranger output
    
    scalefactors_json_file: path to file 'scalefactors_json.json' from spaceranger output, must contain 'spot_diameter_fullres'
    
    downSamplingFactor: some factor that was given (?) to hovernet to reduce image dimensions. Need it here to rescale the coordinates
    and match Space Ranger alignment output files
    
    factor: if 1 then cells under ST spots are considered, multiplies the spot diameter
    
    spot_shape: ['circle', 'square']
        
    space_ranger_output_path: directoriers with two files for each id: /<id>/spatial/scalefactors_json.json and /spatial/tissue_positions_list.csv
    
    hovernet_output_path: directory with files with extract_hov_counts.csv for each id
    
    savepath: directory to save output to
    """
           
    # See data format description at:
    # https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    df_spaceranger_alignment = pd.read_csv(grid_file_path, index_col=0, header=None)
    df_spaceranger_alignment.columns = ['in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
    df_spaceranger_alignment.index.name = 'barcode'

    with open(scalefactors_json_file) as tempfile:
        spot_size_pixels = json.load(tempfile)['spot_diameter_fullres']

    df_hovernet_cells = pd.read_csv(hovernet_data_file_path, index_col=0)
    df_hovernet_cells['cell_type'] = df_hovernet_cells['cell_type'].replace(classes)
    
    func = lambda coords: getSpotID(*coords, df_spaceranger_alignment, spot_diameter_fullres=factor * spot_size_pixels, spot_shape=spot_shape)
    df_hovernet_cells['barcode'] = df_hovernet_cells[['centroid_x', 'centroid_y']].apply(func, axis=1)
    
    if not os.path.exists(savepath):
        os.makedirs(savepath)    
    
    df_hovernet_cells['barcode'].to_csv(savepath + os.path.basename(hovernet_data_file_path) + '.csv')
    
    return df_hovernet_cells


def getSpotID(centroid_y, centroid_x, df_grid_pixels, downSamplingFactor: float = 1., spot_diameter_fullres: float = None, spot_shape: str = 'circle'):
    
    """
    Parameters:
    centroid_x: nucleus location x
    
    centroid_y: nucleus location y
    
    df_grid_pixels: table containing columns 'pxl_row_in_fullres' and 'pxl_col_in_fullres'
    
    downSamplingFactor: factor to scale coordinates
    
    spot_shape: ['circle', 'square']
    
    spot_diameter_fullres: spot diameter in pixels in full resolution WSI
    
    Output:
    Spot identifier or barcode    
    """
    
    centroid_x *= downSamplingFactor
    centroid_y *= downSamplingFactor
    
    se_x = df_grid_pixels['pxl_row_in_fullres']
    se_y = df_grid_pixels['pxl_col_in_fullres']

    if spot_shape=='circle':
        d = np.sqrt((se_x.values - centroid_x)**2 + (se_y.values - centroid_y)**2)
        closest = np.argmin(d)
        if d[closest] <= spot_diameter_fullres/2:
            id =  df_grid_pixels.index[closest]
        else:
            id = None
    
    elif spot_shape=='square':
        in_tile = ((se_x - centroid_x).abs().values <= (spot_diameter_fullres/2)) & \
            ((se_y - centroid_y).abs().values <= (spot_diameter_fullres/2))
        in_tile = np.where(in_tile)[0]
        
        if len(in_tile)>0:
            d = np.sqrt((se_x.values[in_tile] - centroid_x)**2 + (se_y.values[in_tile] - centroid_y)**2)
            closest = np.argmin(d)   
            id =  df_grid_pixels.index[in_tile[closest]]
        else:
            id = None
            
    else:
        raise NotImplementedError
        
    return id


def calculateAggregateValues(df_hovernet_cells, min_cell_type_prob: float = 0.75, savepath: str = None, sname: str = ''):
    
    """
    df_hovernet_cells: pandas.DataFrame with hovernet data assigned to nuclei
    
    min_cell_type_prob: set cell type 'no_label' where probability of assignment is below 'min_cell_type_prob'
    
    savepath: path to save data
    
    sname: file identifier to save data
    """
    
    condition = df_hovernet_cells['cell_type_prob'] < min_cell_type_prob
    df_hovernet_cells.loc[condition, 'cell_type'] = 'no_label'
    
    df_filtered = df_hovernet_cells.loc[~df_hovernet_cells['barcode'].isna()].reset_index().set_index('barcode')

    df_class_densities = pd.get_dummies(df_filtered['cell_type']).groupby(level=0).mean()
    df_class_densities.columns = 'density_' + df_class_densities.columns

    df_class_counts = pd.get_dummies(df_filtered['cell_type']).replace({0: np.nan}).groupby(level=0).count()
    df_class_counts.columns = 'count_' + df_class_counts.columns

    df_total_counts = df_class_counts.sum(axis=1).to_frame()
    df_total_counts.columns = ['total_count']

    df = pd.concat([df_class_densities, df_class_counts, df_total_counts], axis=1)
    
    for q in ['perimeter_length', 'area', 'orientation', 'eccentricity', 'cell_type_prob']:
        # Average per class
        #df_q_classes = []
        #for barcode in df_filtered.index.unique():
        #    df_temp = df_filtered.loc[barcode]
        #    if type(df_temp) is pd.Series:
        #        df_temp = df_temp.to_frame().T
        #    df_q_classes.append(df_temp.set_index('cell_type')[q].astype(float).groupby(level=0).mean().rename(barcode))
        #
        #df_q_classes = pd.concat(df_q_classes, axis=1).T
        df_q_classes = df_filtered.set_index('cell_type', append=True)[q].astype(float).groupby(level=[0, 1]).mean().unstack('cell_type').reindex(df_filtered.index.unique())
        df_q_classes.columns = 'average_%s_' % q + df_q_classes.columns
        df_q_classes.index.name = 'barcode'
        
        # Standard deviation per class
        #df_sdq_classes = []
        #for barcode in df_filtered.index.unique():
        #    df_temp = df_filtered.loc[barcode]
        #    if type(df_temp) is pd.Series:
        #        df_temp = df_temp.to_frame().T
        #    df_sdq_classes.append(df_temp.set_index('cell_type')[q].astype(float).groupby(level=0).std().rename(barcode))
        #
        #df_sdq_classes = pd.concat(df_sdq_classes, axis=1).T
        df_sdq_classes = df_filtered.set_index('cell_type', append=True)[q].astype(float).groupby(level=[0, 1]).std().unstack('cell_type').reindex(df_filtered.index.unique())
        df_sdq_classes.columns = 'std_%s_' % q + df_sdq_classes.columns
        df_sdq_classes.index.name = 'barcode'                 
                
        # Average overall
        df_q = df_filtered[q].groupby(level=0).mean().to_frame()
        df_q.columns = ['average_%s' % q]
        
        # Standard deviation overall
        df_sdq = df_filtered[q].groupby(level=0).std().to_frame()
        df_sdq.columns = ['std_%s' % q]

        df = pd.concat([df, df_q_classes, df_sdq_classes, df_q, df_sdq], axis=1)

    if not savepath is None:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
    
        df.to_csv(savepath + '%s.csv' % sname)

    return df

# Adapted from staintools
def get_stain_matrix(I, angular_percentile=99):
    
    assert I.dtype == np.uint8, "Image should be RGB uint8."
    
    def convert_RGB_to_OD(I):
        mask = (I == 0)
        I[mask] = 1
        return np.maximum(-1 * np.log(I / 255), 1e-6)
    
    OD = convert_RGB_to_OD(I).reshape((-1, 3))

    # Eigenvectors of cov in OD space (orthogonal as cov symmetric)
    _, V = np.linalg.eigh(np.cov(OD, rowvar=False))

    # The two principle eigenvectors
    V = V[:, [2, 1]]

    # Make sure vectors are pointing the right way
    if V[0, 0] < 0: V[:, 0] *= -1
    if V[0, 1] < 0: V[:, 1] *= -1

    # Project on this basis.
    That = np.dot(OD, V)

    # Angular coordinates with repect to the prinicple, orthogonal eigenvectors
    phi = np.arctan2(That[:, 1], That[:, 0])

    # Min and max angles
    minPhi = np.percentile(phi, 100 - angular_percentile)
    maxPhi = np.percentile(phi, angular_percentile)

    # the two principle colors
    v1 = np.dot(V, np.array([np.cos(minPhi), np.sin(minPhi)]))
    v2 = np.dot(V, np.array([np.cos(maxPhi), np.sin(maxPhi)]))

    # Order of H and E.
    # H first row.
    if v1[0] > v2[0]:
        HE = np.array([v1, v2])
    else:
        HE = np.array([v2, v1])
        
    def normalize_matrix_rows(A):
        return A / np.linalg.norm(A, axis=1)[:, None]

    return normalize_matrix_rows(HE)

def calculate_H_E_OD_quantities(img, nuclei, coords, df_stardist, offset='auto', expand_nuclei_distance=15):

    # Calculate cytoplasms mask from nuclei mask
    cytoplasms = expand_labels(nuclei, distance=expand_nuclei_distance)
    cytoplasms[cytoplasms==nuclei] = 0

    HE = get_stain_matrix(img)
    print('Hematoxylin:\t', HE[0])
    print('Eosin:\t\t', HE[1])
    H, E = np.moveaxis(np.dot(img, HE.T), 2, 0)

    H /= H.max()
    E /= E.max()
    
    if offset == 'auto':
        N = 5000
        se = df_stardist['area'].sample(N, replace=False) if df_stardist.shape[0] > N else df_stardist['area'].copy()
        offset = int(se.apply(np.sqrt).quantile(0.99) * 2)
    print('Offset:', offset)
    
    print(coords)
    
    d = dict()
    ids = np.unique(np.ravel(nuclei))
    print(len(ids))
    for id in ids[:]:
        if id != 0:
            sid = str(id - 1)
            d[sid] = dict()

            center_raw = df_stardist.loc[sid][['centroid_x', 'centroid_y']].astype(int).values
            center = center_raw[0] - coords[1], center_raw[1] - coords[0]
            x_min, x_max = max(0, center[1] - offset), min(H.shape[0], center[1] + offset)
            y_min, y_max = max(0, center[0] - offset), min(H.shape[1], center[0] + offset)
            
            wh_cy = np.where(cytoplasms[x_min: x_max, y_min: y_max]==id)
            wh_cy = wh_cy[0]  + x_min, wh_cy[1]  + y_min
            wh_nuc = np.where(nuclei[x_min: x_max, y_min: y_max]==id)
            wh_nuc = wh_nuc[0]  + x_min, wh_nuc[1]  + y_min

            if (len(wh_cy[0]) > 0) and (len(wh_nuc[0]) > 0):
                d[sid].update({'x': center_raw[0],
                              'y': center_raw[1],
                              'cy_hemat_mean': H[wh_cy].mean(),
                              'cy_eosin_mean': E[wh_cy].mean(),
                              'nuc_hemat_mean': H[wh_nuc].mean(),
                              'nuc_eosin_mean': E[wh_nuc].mean(),
                              'cy_eosin_max': E[wh_cy].max(),
                              'cy_eosin_min': E[wh_cy].min()})

    return pd.DataFrame(d).T

