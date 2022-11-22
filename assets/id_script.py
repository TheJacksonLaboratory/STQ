
ids_WM4237_st = {'WM4237_TE_S1_ST': 'SC2200092_WM4237-1-1337',
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

print('sample,image,grid')
for id in sorted(ids_WM4237_st.keys()):
    print('%s,/projects/chuang-lab/rubinstein/ST_melanomaPDX/data/%s/img/%s.tiff,/projects/chuang-lab/rubinstein/ST_melanomaPDX/data/%s/spaceranger/spatial/' \
          % (id, ids_WM4237_st[id], ids_WM4237_st[id].split('_')[0], ids_WM4237_st[id]))
          
 



         
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
                 'WM4007_T4_S2_ST': 'SC2200811_WM4007-T4-6',
                }

print('sample,image,grid')
for id in sorted(ids_WM4007_st.keys()):
    print('%s,/projects/chuang-lab/rubinstein/ST_melanomaPDX/data/WM4007/%s/img/%s.tiff,/projects/chuang-lab/rubinstein/ST_melanomaPDX/data/WM4007/%s/spaceranger/spatial/' \
          % (id, ids_WM4007_st[id], ids_WM4007_st[id].split('_')[0], ids_WM4007_st[id]))    )    







ids_WM4007_add ={'WM4007_TC_S1_AD': 'CTRL_#439_6-28-22__20221020_081352_cut_level_0_oid_0',
                 'WM4007_TC_S2_AD': 'CTRL_#439_8-18-22_20221020_080124_cut_level_0_oid_0',
                 'WM4007_T0_S1_AD': 'TO_#484_6_28_22__20221020_081152_cut_level_0_oid_1',
                 'WM4007_T0_S2_AD': 'TO_484_8-18-22_20221020_080241_cut_level_0_oid_0',
                 'WM4007_T1_S1_AD': 'T1_#426_6-28-22__20221020_081232_cut_level_0_oid_0',
                 'WM4007_T1_S2_AD': 'T1_426_8-18-22_20221020_080340_cut_level_0_oid_0',
                 'WM4007_T2_S1_AD': 'T2_#437_6-28-22__20221020_081518_cut_level_0_oid_0',
                 'WM4007_T2_S2_AD': 'T2_437_8-18-22_20221020_080440_cut_level_0_oid_0',
                 'WM4007_T3_S1_AD': 'T3_#17_6-28-22__20221020_081603_cut_level_0_oid_0',
                 'WM4007_T3_S2_AD': 'T3_17_8-18-22_20221020_080910_cut_level_0_oid_0',
                 'WM4007_T4_S1_AD': 'T4_#006_6_28_22__20221020_081656_cut_level_0_oid_0',
                 'WM4007_T4_S2_AD': 'T4_6_8-18-22_20221020_081100_cut_level_0_oid_0'}

print('sample,image,grid')
for id in sorted(ids_WM4007_add.keys()):
    print('%s,/projects/chuang-lab/USERS/domans/melanoma_PDX_ST/additionalHandE/tiff/WM4007/%s.tiff,' % (id, ids_WM4007_add[id]))




 
ids_WM4237_add={"WM4237_TE_S1_AD": "WM4237 CTRL_200__20220520_084527_cut_level_0_oid_1",
                "WM4237_TE_S2_AD": "WM4237 CTRL_#200__20220520_084919_cut_level_0_oid_1",
                "WM4237_T0_S1_AD": "WM4237 #197_TO__20220520_082852_cut_level_0_oid_2",
                "WM4237_T0_S2_AD": "WM4237 197_T0__20220520_091543_cut_level_0_oid_0",
                "WM4237_TC_S1_AD": "WM4237 CTRL_207_DAY_42__20220520_084211_cut_level_0_oid_0",
                "WM4237_TC_S2_AD": "WM4237 CTRL_207__20220520_085512_cut_level_0_oid_0",
                "WM4237_T1_S1_AD": "WM4237 #229_T1_Day_14__20220520_083037_cut_level_0_oid_1",
                "WM4237_T1_S2_AD": "WM4237 #229_T1__20220520_091746_cut_level_0_oid_0",
                "WM4237_T2_S1_AD": "WM4237 #211_T2_Day_42__20220520_083124_cut_level_0_oid_0",
                "WM4237_T2_S2_AD": "WM4237 #211_T2__20220520_085201_cut_level_0_oid_1",
                "WM4237_T3_S1_AD": "WM4237 #221_T3_Day_70__20220520_083633_cut_level_0_oid_0",
                "WM4237_T3_S2_AD": "WM4237 #221_T3__20220520_085259_cut_level_0_oid_0",
                "WM4237_T4_S1_AD": "WM4237 #214_T4_Day_94__20220520_084017_cut_level_0_oid_0",
                "WM4237_T4_S2_AD": "WM4237 #214_T4__20220520_085414_cut_level_0_oid_0"}

{k:ids_WM4237_add[k] for k in sorted(ids_WM4237_add.keys())}

print('sample,image,grid')
for id in sorted(ids_WM4237_add.keys()):
    print('%s,/projects/chuang-lab/USERS/domans/melanoma_PDX_ST/additionalHandE/tiff/WM4007/%s.tiff,' % (id, ids_WM4237_add[id]))
    
