import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

def prepInfoJSON(tdatapath, model, fname, eps=50, min_samples=500, selobj=None, f=1.5, efx=0.35, efy=0.35, manshift=[0,0,0,0]):

    img_RGB_high_res = plt.imread(tdatapath + fname)

    img_RGB_high_res  = img_RGB_high_res[::3, ::3, :]

    v = img_RGB_high_res[:, :, :3].copy()
    v = v.mean(axis=2)
    vc = v.copy()
    v[vc<100] = 255
    v[vc>200] = 255
    v[(vc>=100) & (vc<=200)] = 0
    v = pd.DataFrame(v.T).replace({255: np.nan}).stack().dropna().index.to_frame().values

    np.random.seed(0)
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(v)
    labels = db.labels_
    vco = pd.Series(labels).value_counts()
    
    fig, axs = plt.subplots(1, 2, figsize=(f*10, f*4))

    axs[0].imshow(img_RGB_high_res[:, :, :3], origin='lower')
    axs[0].set_xticks([])
    axs[0].set_xticklabels([])
    axs[0].set_yticks([])
    axs[0].set_yticklabels([])
    
    axs[1].scatter(v.T[0], v.T[1], c=labels, s=1)
    for l in vco.index:
        if l!=-1:
            axs[1].text(np.median(v.T[0][np.where(labels==l)[0]]),
                     np.median(v.T[1][np.where(labels==l)[0]]),
                     l, va='center', ha='center', fontsize=20)

    id = fname[:-len('_level_6.tiff')]

    objects = sorted(list(set(vco.to_dict().keys()).difference([-1])))
    res = dict()
    for l in objects:
        res[l] = dict()
        res[l]['model'] = model
        res[l]['id'] = id
        for i in range(2):
            res[l][i] = dict()
            se_temp = pd.Series(v.T[i][np.where(labels==l)[0]])
            a, b = se_temp.quantile(0.05), se_temp.quantile(0.95)
            sp = b - a

            ef = efx if i==0 else efy
            sa, sb = (manshift[0], manshift[1]) if i==0 else (manshift[2], manshift[3])
            a -= (ef+sa)*sp
            b += (ef+sb)*sp

            d = img_RGB_high_res.shape[1 if i==0 else 0]
            a, b = a/d, b/d - a/d
            a = max(a, 0)
            if a + b > 1:
                b = 1 - a
            res[l][i]['location'] = a
            res[l][i]['size'] = b
    
    axs[0].text(axs[0].get_xlim()[0], axs[0].get_ylim()[0], id, va='bottom', ha='left')
    
    axs[1].set_xticks([])
    axs[1].set_xticklabels([])
    axs[1].set_yticks([])
    axs[1].set_yticklabels([])
        
    # Save to JSON
    if not selobj is None:
        with open(tdatapath + '%s_%s_%s.json' % (model, id, selobj), 'w') as outfile:
            outfile.write(json.dumps(res[selobj]))
        axs[1].text(axs[1].get_xlim()[0], axs[1].get_ylim()[0], 'Selected object: %s' % selobj, va='bottom', ha='left')
        
        x1 =  res[selobj][0]['location']*img_RGB_high_res.shape[1]
        x2 =  x1 + res[selobj][0]['size']*img_RGB_high_res.shape[1]
        y1 =  res[selobj][1]['location']*img_RGB_high_res.shape[0]
        y2 =  y1 + res[selobj][1]['size']*img_RGB_high_res.shape[0]
        axs[1].plot([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], linewidth=0.5, c='crimson')

        print(res[selobj][0]['size'] * res[selobj][1]['size'])
           
    axs[1].set_aspect('equal')
    plt.show()
    
    return

