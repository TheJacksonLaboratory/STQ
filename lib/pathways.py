
import numpy as np
import pandas as pd

genesDicts = dict()
loadPathwayRGD = lambda name: np.sort(pd.read_csv(name, delimiter='\t', index_col=0)['Symbol'].values.tolist())
loadPathway = lambda name: np.sort(pd.read_csv(name, delimiter='\t', index_col=0).loc['MAPPED_SYMBOLS'].values).tolist()[0].split(',')

filesPath = 'c:/Projects/A_ST/'

genesDicts.update({'PKA': loadPathwayRGD(filesPath + 'annotation PKA.tab')})

genesDicts.update({'MAPK': loadPathway(filesPath + 'KEGG_MAPK_SIGNALING_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'AKT': loadPathway(filesPath + 'REACTOME_PI3K_AKT_SIGNALING_IN_CANCER.v7.5.1.tsv')})
genesDicts.update({'TRKA': loadPathway(filesPath + 'BIOCARTA_TRKA_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'MELANOMA': loadPathway(filesPath + 'KEGG_MELANOMA.v7.5.1.tsv')})
genesDicts.update({'CELL_CYCLE': loadPathway(filesPath + 'KEGG_CELL_CYCLE.v7.5.1.tsv')})
genesDicts.update({'MELANOGENESIS': loadPathway(filesPath + 'KEGG_MELANOGENESIS.v7.5.1.tsv')})
genesDicts.update({'PI3K_AKT_MTOR': loadPathway(filesPath + 'HALLMARK_PI3K_AKT_MTOR_SIGNALING.v7.5.1.tsv')})
genesDicts.update({'REACTOME_MICRORNA': loadPathway(filesPath + 'REACTOME_MICRORNA_MIRNA_BIOGENESIS.v7.5.1.tsv')})
genesDicts.update({'VEGFAVEGFR2': loadPathway(filesPath + 'WP_VEGFAVEGFR2_SIGNALING_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'AP1': loadPathway(filesPath + 'PID_AP1_PATHWAY.v7.5.1.tsv')})


genesDicts.update({'KEGG_FOCAL_ADHESION': loadPathway(filesPath + 'KEGG_FOCAL_ADHESION.v2022.1.Hs.tsv')})
genesDicts.update({'GOBP_CELL_MOTILITY': loadPathway(filesPath + 'GOBP_CELL_MOTILITY.v2022.1.Hs.tsv')})


genesDicts.update({'ERK': loadPathway(filesPath + 'BIOCARTA_ERK_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'MAPK1_ERK2': loadPathway(filesPath + 'REACTOME_MAPK1_ERK2_ACTIVATION.v7.5.1.tsv')})
genesDicts.update({'MAPK3_ERK1': loadPathway(filesPath + 'REACTOME_MAPK3_ERK1_ACTIVATION.v7.5.1.tsv')})
genesDicts.update({'MTOR': loadPathway(filesPath + 'KEGG_MTOR_SIGNALING_PATHWAY.v7.5.1.tsv')})

genesDicts.update({'JAK_STAT': loadPathway(filesPath + 'JAK_STAT_CASCADE.v7.5.1.tsv')})
genesDicts.update({'KEGG_JAK_STAT': loadPathway(filesPath + 'KEGG_JAK_STAT_SIGNALING_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'WP_HIPPO': loadPathway(filesPath + 'WP_HIPPO_SIGNALING_REGULATION_PATHWAYS.v7.5.1.tsv')})
genesDicts.update({'HIPPO': loadPathway(filesPath + 'REACTOME_SIGNALING_BY_HIPPO.v7.5.1.tsv')})
genesDicts.update({'NRF2': loadPathway(filesPath + 'WP_NRF2_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'FAO': loadPathway(filesPath + 'FATTY_ACID_OXIDATION.v7.5.1.tsv')})

genesDicts.update({'KEGG_FATTY_ACID': loadPathway(filesPath + 'KEGG_FATTY_ACID_METABOLISM.v2022.1.Hs.tsv')})

genesDicts.update({'Reactome-FAO': loadPathway(filesPath + 'REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION.v7.5.1.tsv')})
genesDicts.update({'OXPHOS': loadPathway(filesPath + 'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA.v7.5.1.tsv')})
genesDicts.update({'FAK': loadPathway(filesPath + 'PID_FAK_PATHWAY.v7.5.1.tsv')})
genesDicts.update({'HEAT': loadPathway(filesPath + 'REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS.v7.5.1.tsv')})
genesDicts.update({'NEURAL_CREST': loadPathway(filesPath + 'WP_NEURAL_CREST_DIFFERENTIATION.v7.5.1.tsv')})
genesDicts.update({'KEGG_LYSOSOME': loadPathway(filesPath + 'KEGG_LYSOSOME.v2022.1.Hs.tsv')})
genesDicts.update({'HALLMARK_APOPTOSIS': loadPathway(filesPath + 'HALLMARK_APOPTOSIS.v2022.1.Hs.tsv')})


genesDicts.update({'G1_S_CELL_CYCLE': loadPathway(filesPath + 'FISCHER_G1_S_CELL_CYCLE.v7.5.1.tsv')})
genesDicts.update({'G2_M_CELL_CYCLE': loadPathway(filesPath + 'FISCHER_G2_M_CELL_CYCLE.v7.5.1.tsv')})


genesDicts.update({'GOBP_MELANOCYTE_PROLIFERATION': loadPathway(filesPath + 'GOBP_MELANOCYTE_PROLIFERATION.v2022.1.Hs.tsv')})
genesDicts.update({'GOBP_MESEN_PROLIFERATION': loadPathway(filesPath + 'GOBP_MESENCHYMAL_CELL_PROLIFERATION.v2022.1.Hs.tsv')})
genesDicts.update({'REACTOME_VEGFR2_PROLIFERATION': loadPathway(filesPath + 'REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION.v2022.1.Hs.tsv')})
genesDicts.update({'WP_GLYCOLYSIS_AND_GLUCONEOGENESIS': loadPathway(filesPath + 'WP_GLYCOLYSIS_AND_GLUCONEOGENESIS.v2022.1.Hs.tsv')})
genesDicts.update({'WP_ONECARBON_METABOLISM': loadPathway(filesPath + 'WP_ONECARBON_METABOLISM.v2022.1.Hs.tsv')})
genesDicts.update({'HALLMARK_EMT': loadPathway(filesPath + 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2022.1.Hs.tsv')})


genesDicts.update({'BIOCARTA_P53HYPOXIA_PATHWAY': loadPathway(filesPath + 'BIOCARTA_P53HYPOXIA_PATHWAY.v2022.1.Hs.tsv')})
genesDicts.update({'RESPONSE_TO_HYPOXIA': loadPathway(filesPath + 'RESPONSE_TO_HYPOXIA.v2022.1.Hs.tsv')})
genesDicts.update({'PID_HIF1A_PATHWAY': loadPathway(filesPath + 'PID_HIF1A_PATHWAY.v2022.1.Hs.tsv')})
genesDicts.update({'PID_HIF1_TFPATHWAY': loadPathway(filesPath + 'PID_HIF1_TFPATHWAY.v2022.1.Hs.tsv')})


genesDicts.update({'OREXIN_SIG': 'CCN2,DUSP4,EGR1,EGR3,FOSL2,NR4A1,JUNB,MMP13,RARA,SPP1,NR4A3,BHLHE40,CXXC5,LIF,MCAM,NGFR,FOSL1,SPHK1,RAMP1,CEMIP2,TNFRSF12A,MET,SEMA4C,CCN3,TGFBI'.split(',')})
genesDicts.update({'NUCLEAR_SIG': 'DUSP4,DUSP6,EGR1,EGR3,JUNB,MEF2D,VGF,FOSL1,TRIB1,BRAF,COL6A3,FGFR1,MET,SPHK1,SH2B3,CCND1,FOSL2,LIF,NR4A1,TFAP2A,B2M,NGFR,RARA,RCAN1'.split(',')})
genesDicts.update({'REG_MAPK_SIG': 'BRAF,CCN2,DUSP4,DUSP5,DUSP6,FGFR1,LIF,MAP3K11,TRAF1,SPHK1,SH2B3,DUSP10,SEMA4C,SPRY4,DUSP15,NR4A1,MET,NGFR,DNAJB1,MPRIP,AGK,PTPRF,GOLGA4,JUNB,PSMA2,TRIB1,RCAN1,RABGEF1,PPP1R15A,PFKFB4,TNS3'.split(',')})



SenMayo = 'ACVR1B,ANG,ANGPT1,ANGPTL4,AREG,AXL,BEX3,BMP2,BMP6,C3,CCL1,CCL13,CCL16,CCL2,CCL20,CCL24,CCL26,CCL3,CCL3L1,CCL4,CCL5,CCL7,CCL8,CD55,CD9,CSF1,CSF2,CSF2RB,CST4,CTNNB1,CTSB,CXCL1,CXCL10,CXCL12,CXCL16,CXCL2,CXCL3,CXCL8,CXCR2,DKK1,EDN1,EGF,EGFR,EREG,ESM1,ETS2,FAS,FGF1,FGF2,FGF7,GDF15,GEM,GMFG,HGF,HMGB1,ICAM1,ICAM3,IGF1,IGFBP1,IGFBP2,IGFBP3,IGFBP4,IGFBP5,IGFBP6,IGFBP7,IL10,IL13,IL15,IL18,IL1A,IL1B,IL2,IL32,IL6,IL6ST,IL7,INHA,IQGAP2,ITGA2,ITPKA,JUN,KITLG,LCP1,MIF,MMP1,MMP10,MMP12,MMP13,MMP14,MMP2,MMP3,MMP9,NAP1L4,NRG1,PAPPA,PECAM1,PGF,PIGF,PLAT,PLAU,PLAUR,PTBP1,PTGER2,PTGES,RPS6KA5,SCAMP4,SELPLG,SEMA3F,SERPINB4,SERPINE1,SERPINE2,SPP1,SPX,TIMP2,TNF,TNFRSF10C,TNFRSF11B,TNFRSF1A,TNFRSF1B,TUBGCP2,VEGFA,VEGFC,VGF,WNT16,WNT2'
genesDicts.update({'SenMayo': SenMayo.split(',')})

genesDicts.update({'hypoxia26': 'ALDOA,ANGPTL4,ANLN,BNC1,C20orf20,CA9,CDKN3,COL4A6,DCBLD1,ENO1,FAM83B,FOSL1,GNAI1,HIG2,KCTD11,KRT17,LDHA,MPRS17,P4HA1,PGAM1,PGK1,SDC1,SLC16A1,SLC2A1,TPI1,VEGFA'.split(',')})

# KEGG_JAK_STAT!!!
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3274382/
genesDicts.update({'ISG1': ['ADAR','APOBEC3','BST2','C6orf150','CD74','DDIT4','DDX58','DDX60','EIF2AK2','GBP1','HPSE','IFI44L','IFI6','IFIH1','IFIT1','IFITM1',
       'IRF1','IRF7','ISG15','ISG20','MAP3K14','MOV10','MS4A4A','MX1','MX2','NAMPT','NT5C3','OAS1','OASL','P2RY6','PHF15','PML','RSAD2','RTP4',
       'SLC15A3','SLC25A28','SSBP3','TREX1','TRIM5','TRIM25','SUN2','ZC3HAV1','MB21D1','PKR','GBP2','G1P3','MDA5','IFIT2','IFITM2','IFIT3',
       'IFITM5','IFIT5','NIK','PBEF1','OAS2','OAS3','TRIM19','AGS1','UNC84B','ZAP']})

temp = 'Casp1, Casp4, Ccl4, Cd274, Cmpk2, Cxcl10, Gbp2, Gbp6, Herc6, Ifi204,Ccl13, Ccl4, Cxcl10, Ddx58, Ifih1, Ifit1b, Irf7, Rsad2, Zbp1,Casp1, Ccl4, Cxcl10, Gbp2, Ifi204, Ifi35, Ifit1b, Irf1, Psmb8, Psmb9,Ccl13, Cxcl10, Gbp2, Ifi204, Ifit3, Irf1, Isg20, Lgals3bp, Ly6a, Mov10,Ifit1b, Ifit3, Irf1, Irf7, Irf9,Irf1, Stat1, Stat2, Tgtp1,Ccl4, Cxcl10, Ifit1b, Irf1, Nmi, Stat1,Cxcl10, Ifih1, Ifit1b, Isg20, Trim30,Cxcl10, Irf1, Irf7,Ifit1b, Zbp1,Casp1, Gbp2, Il18bp,'
temp = temp.split(',')
temp = [gene.strip().upper() for gene in temp]
genesDicts.update({'ISG2': temp})

genesDicts.update({'ISG3': ['STAT1','IRF3','STAT3','STAT6','STAT2','IRF9','NFKB1, NFKB2, REL, RELA, RELB','STAT4','IRF2','IRF7','IRF1']})

# GABPA is NRF2