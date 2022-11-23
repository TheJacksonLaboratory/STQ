
ids_WM4237_st ={'WM4237_TE_S1_ST': 'SC2200092_WM4237-1-1337',
                'WM4237_TE_S2_ST': 'SC2200093_WM4237-1-1337',
                'WM4237_TE_S3_ST': 'SC2200094_WM4237-1-1337',
                'WM4237_TE_S4_ST': 'SC2200095_WM4237-1-1337',
                'WM4237_T0_S1_ST': 'SC2200319_WM4237-T0-Day0-197',
                'WM4237_T0_S2_ST': 'SC2200320_WM4237-T0-Day0-197',
                'WM4237_TC_S1_ST': 'SC2200321_WM4237-CTRL-Day42-207',
                'WM4237_TC_S2_ST': 'SC2200322_WM4237-CTRL-Day42-207',
                'WM4237_T1_S1_ST': 'SC2200323_WM4237-T1-Day14-229',
                'WM4237_T1_S2_ST': 'SC2200324_WM4237-T1-Day14-229',
                'WM4237_T2_S1_ST': 'SC2200325_WM4237-T2-Day42-211', 
                'WM4237_T2_S2_ST': 'SC2200326_WM4237-T2-Day42-211',
                'WM4237_T3_S1_ST': 'SC2200327_WM4237-T3-Day70-221',
                'WM4237_T3_S2_ST': 'SC2200328_WM4237-T3-Day70-221',
                'WM4237_T4_S1_ST': 'SC2200329_WM4237-T4-Day94-214',
                'WM4237_T4_S2_ST': 'SC2200330_WM4237-T4-Day94-214'}

data = ['sample,fastq,image\n']
basepath = '/projects/chuang-lab/rubinstein/ST_melanomaPDX/data'
tids = ids_WM4237_st.copy()
for i, id in enumerate(sorted(tids.keys())):
    s = '%s,%s/%s/fastq/,%s/%s/img/%s.tiff' \
          % (id, basepath, tids[id], basepath, tids[id], tids[id].split('_')[0])
    if i < len(tids) - 1:
        s += '\n'
    data.append(s)
with open('samplesheet_WM4237.csv', 'w') as tempfile:
    tempfile.writelines(data)
    

                
ids_WM4007_st = {'WM4007_TC_S1_ST': 'SC2200800_WM4007-Ctrl-439',
                 'WM4007_TC_S2_ST': 'SC2200801_WM4007-Ctrl-439',                
                 'WM4007_T0_S1_ST': 'SC2200802_WM4007-T0-484',
                 'WM4007_T0_S2_ST': 'SC2200803_WM4007-T0-484',     
                 'WM4007_T1_S1_ST': 'SC2200804_WM4007-T1-426',
                 'WM4007_T1_S2_ST': 'SC2200805_WM4007-T1-426',               
                 'WM4007_T2_S1_ST': 'SC2200806_WM4007-T2-437',
                 'WM4007_T2_S2_ST': 'SC2200807_WM4007-T2-437',                 
                 'WM4007_T3_S1_ST': 'SC2200808_WM4007-T3-17',
                 'WM4007_T3_S2_ST': 'SC2200809_WM4007-T3-17',                
                 'WM4007_T4_S1_ST': 'SC2200810_WM4007-T4-6',
                 'WM4007_T4_S2_ST': 'SC2200811_WM4007-T4-6'}

data = ['sample,fastq,image\n']
basepath = '/projects/chuang-lab/rubinstein/ST_melanomaPDX/data/WM4007'
tids = ids_WM4007_st.copy()
for i, id in enumerate(sorted(tids.keys())):
    s = '%s,%s/%s/fastq/,%s/%s/img/%s.tiff' \
          % (id, basepath, tids[id], basepath, tids[id], tids[id].split('_')[0])
    if i < len(tids) - 1:
        s += '\n'
    data.append(s)
with open('samplesheet_WM4007.csv', 'w') as tempfile:
    tempfile.writelines(data)
    