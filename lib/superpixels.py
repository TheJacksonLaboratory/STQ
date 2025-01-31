
import json
import gzip
import numpy as np
import cv2
import matplotlib.pyplot as plt

def plot_spx_contours(all_contours, figsize=(15, 15)):
    fig, ax = plt.subplots(figsize=figsize)
    for c in all_contours.keys():
        xp = []
        yp = []
        for sub_contour in all_contours[c]:
            x = np.array(sub_contour).T[0].tolist()
            y = np.array(sub_contour).T[1].tolist()
            xp += x + [x[0], None]
            yp += y + [y[0], None]
        ax.plot(xp, yp, '-o', ms=0, lw=2, label=c)
    if len(all_contours.keys()) < 10:
        plt.legend()
    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()
    return

def get_countours_from_mask(superpixelation):

    '''
    get_countours_from_mask(np.array([[0, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
                                      [0, 0, 2, 2, 2, 2, 2, 2, 2, 8, 8, 8],
                                      [0, 0, 2, 2, 5, 5, 2, 2, 2, 8, 8, 8],
                                      [0, 0, 2, 2, 5, 5, 2, 2, 2, 8, 8, 8],
                                      [0, 0, 2, 2, 2, 2, 2, 2, 2, 8, 8, 8],
                                      [0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 8, 8],
                                      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3],
                                      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3],
                                      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3],
                                      [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3]]))
    '''

    raw_contours = cv2.findContours(superpixelation + 1, cv2.RETR_FLOODFILL, cv2.CHAIN_APPROX_SIMPLE)
    all_contours = dict()
    for c in range(len(raw_contours[0])):
        cid = str(superpixelation[raw_contours[0][c][0][0][1], raw_contours[0][c][0][0][0]])
        if not cid in all_contours.keys():
            all_contours.update({cid: []})
        new_contour = raw_contours[0][c][:, 0, :].tolist()
        all_contours.update({cid: all_contours[cid] + [new_contour]})
    
    print('Created %s contours:' % len(all_contours))
    
    return all_contours

def save_contours(all_contours, filename='contours.json.gz'):

    with gzip.GzipFile(filename, 'w') as tempfile:
        tempfile.write(json.dumps(all_contours).encode('utf-8'))
    
    return