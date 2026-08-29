import numpy as np
import cv2
from scipy.ndimage import find_objects

def get_contours(mask):
    """Extract contours from uint32 segmentation mask.
    """
    contours = {}
    centroids = {}
    slices = find_objects(mask)
    
    for cell_id, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        # Extract cell region and create binary mask
        region = mask[slc] == cell_id
        # Find contours (returns list, take the largest)
        cnts, _ = cv2.findContours(region.astype(np.uint8), 
                                    cv2.RETR_EXTERNAL, 
                                    cv2.CHAIN_APPROX_NONE)
        if cnts:
            # Offset contour coordinates to original image space
            cnt = max(cnts, key=cv2.contourArea)
            cnt[:, 0, 0] += slc[1].start
            cnt[:, 0, 1] += slc[0].start
            contours[cell_id] = cnt[:, 0, :]

            p = cnt[:, 0, :]
            y, x = p[:, 0], p[:, 1]
            x2, y2 = np.roll(x, -1), np.roll(y, -1)
            cross = x * y2 - x2 * y
            A = cross.sum()
            centroid = np.array([
                (y + y2).dot(cross) / (3 * A),
                (x + x2).dot(cross) / (3 * A)])
            if centroid[0]!= centroid[0] or centroid[1] != centroid[1]:
                centroid = cnt.mean(axis=0)[0]
            
            centroids[cell_id] = centroid

    return contours, centroids
