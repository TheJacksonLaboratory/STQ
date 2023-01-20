import os
import argparse
import pandas as pd
import numpy as np
import tensorflow as tf
from tqdm import tqdm
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute Inception V3 features on tiles')
    parser.add_argument('--input-path', dest='input_path', action='store',
                        required=True,
                        help="""The path to the tiles from a whole slide image (WSI) in a TIF format.""")
    parser.add_argument('--output-path', dest='output_path', action='store',
                        required=True,
                        help="""Name of CSV file in which to store the feature matrix (rows are tiles, cols are features). 
                                The file will be compressed if it is named *.gz""")
    args = parser.parse_args()
     
    output_path = args.output_path
    input_path = args.input_path
    
    # Make tile names list
    fnames = [fname for fname in os.listdir(input_path) if fname[-len('.tif'):]=='.tif']
    num_images = len(fnames)
    print("Number of images:", num_images)
    
    # Assuming that all tiles have the same shape
    # Read the first tile and create a model
    tile = plt.imread(input_path + fnames[0])
   
    base_model = tf.keras.applications.inception_v3.InceptionV3(include_top=False, weights='imagenet', input_shape=tile.shape)
    xi = base_model.output
    xi = tf.keras.layers.GlobalAveragePooling2D(data_format=None)(xi)
    model = tf.keras.models.Model(inputs=base_model.input, outputs=xi)
  
    batch_size = int(10**8 / (tile.shape[0] * tile.shape[1]))
    num_batches = int(np.ceil(num_images / batch_size))
    
    print('Reading and pocessing tiles:', num_images)
    print('Batch size:', batch_size)
    print('Number of batches:', num_batches)
    
    features = []
    for ibatch in tqdm(range(num_batches)):
        images = []
        for indx in range(batch_size):
            try:
                images.append(plt.imread(input_path + fnames[indx + ibatch*batch_size]))
            except:
                pass

        features.append(model.predict(tf.keras.applications.inception_v3.preprocess_input(np.stack(images)), verbose=0))
    features = np.vstack(features)
    
    df_features = pd.DataFrame(data=features, index=[fname[:-len('.tif')] for fname in fnames], columns=['feat' + str(i) for i in range(features.shape[1])])
    print(df_features)
    
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    
    df_features.to_csv(output_path)
    print('Successfully wrote:' + output_path)
    
exit(0)
