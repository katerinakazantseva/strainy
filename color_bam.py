import csv
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt




edges=['edge_1', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_108', 'edge_109', 'edge_110', 'edge_114', 'edge_115', 'edge_116', 'edge_117', 'edge_118', 'edge_120', 'edge_121', 'edge_122', 'edge_123', 'edge_131', 'edge_138', 'edge_139', 'edge_140', 'edge_141', 'edge_142', 'edge_143', 'edge_144', 'edge_145', 'edge_146', 'edge_149', 'edge_194', 'edge_195', 'edge_196', 'edge_197', 'edge_203', 'edge_219', 'edge_22', 'edge_220', 'edge_221', 'edge_231', 'edge_232', 'edge_233', 'edge_234', 'edge_235', 'edge_236', 'edge_237', 'edge_238', 'edge_239', 'edge_240', 'edge_241', 'edge_242', 'edge_249', 'edge_250', 'edge_253', 'edge_254', 'edge_255', 'edge_256', 'edge_257', 'edge_258', 'edge_259', 'edge_260', 'edge_266', 'edge_267', 'edge_269', 'edge_27', 'edge_270', 'edge_271', 'edge_272', 'edge_273', 'edge_275', 'edge_276', 'edge_278', 'edge_279', 'edge_28', 'edge_280', 'edge_281', 'edge_283', 'edge_284', 'edge_285', 'edge_29', 'edge_291', 'edge_294', 'edge_295', 'edge_296', 'edge_3', 'edge_30', 'edge_300', 'edge_301', 'edge_302', 'edge_303', 'edge_304', 'edge_305', 'edge_306', 'edge_307', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_57', 'edge_64', 'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_98', 'edge_99']
edges=['edge_1', 'edge_10', 'edge_100', 'edge_101', 'edge_102', 'edge_103', 'edge_104', 'edge_105', 'edge_106', 'edge_107', 'edge_108', 'edge_109', 'edge_11', 'edge_110', 'edge_111', 'edge_112', 'edge_113', 'edge_114', 'edge_115', 'edge_116', 'edge_12', 'edge_13', 'edge_14', 'edge_15', 'edge_16', 'edge_17', 'edge_18', 'edge_19', 'edge_2', 'edge_20', 'edge_21', 'edge_22', 'edge_23', 'edge_24', 'edge_25', 'edge_26', 'edge_27', 'edge_28', 'edge_29', 'edge_3', 'edge_30', 'edge_31', 'edge_32', 'edge_33', 'edge_34', 'edge_35', 'edge_36', 'edge_37', 'edge_38', 'edge_39', 'edge_4', 'edge_40', 'edge_41', 'edge_42', 'edge_43', 'edge_44', 'edge_45', 'edge_46', 'edge_47', 'edge_48', 'edge_49', 'edge_5', 'edge_50', 'edge_51', 'edge_52', 'edge_53', 'edge_54', 'edge_55', 'edge_56', 'edge_57', 'edge_58', 'edge_59', 'edge_6', 'edge_60', 'edge_61', 'edge_62', 'edge_63', 'edge_64', 'edge_65', 'edge_66', 'edge_67', 'edge_68', 'edge_69', 'edge_7', 'edge_70', 'edge_71', 'edge_72', 'edge_73', 'edge_74', 'edge_75', 'edge_76', 'edge_77', 'edge_78', 'edge_79', 'edge_8', 'edge_80', 'edge_81', 'edge_82', 'edge_83', 'edge_84', 'edge_85', 'edge_86', 'edge_87', 'edge_88', 'edge_89', 'edge_9', 'edge_90', 'edge_91', 'edge_92', 'edge_93', 'edge_94', 'edge_95', 'edge_96', 'edge_97', 'edge_98', 'edge_99']

edges=['edge_100']
I=1000
AF=0.1





#bam ="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/sim.3ecoliY.bam"
#bam="/Users/ekaterina.kazantseva/MT/test_data/SRR13128014_labeled.bam"
bam="/Users/ekaterina.kazantseva/MT/5strain_hifi_sim/3/noalt/sim.3ecoliv2.bam"

def write_bam(bam, edge, I, AF):
    print("color bam")
    #infile = pysam.AlignmentFile(bam, "rb")
    #outfile = pysam.AlignmentFile("output/bam/coloredBAM_%s_%s_%s.bam" % (edge, I, AF), "wb", template=infile)
    cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF),keep_default_na=False)

    iter = infile.fetch(edge,until_eof=True)
    cmap = plt.get_cmap('viridis')
    cl.loc[cl['Cluster'] == 'NA', 'Cluster'] = 0
    clusters=sorted(set(cl['Cluster'].astype(int)))
    cmap = cmap(np.linspace(0, 1, len(clusters)))
    colors={}
    i=0
    colors[0] = "#505050"
    try:
        clusters.remove('0')
    except: KeyError



    for cluster in clusters:
        colors[cluster]=mt.colors.to_hex(cmap[i])
        i=i+1
    cl_dict = dict(zip(cl.ReadName, cl.Cluster))

    print(colors)
    #print(cl_dict)
    for read in iter:

        try:

            clN=int(cl_dict[str(read).split()[0]])

            tag=colors[clN]
            #print(clN)
            #print(tag)
            #print(str(read).split()[0]+" "+str(clN)+" "+str(tag))
            read.set_tag('YC', tag, replace=False)


            outfile.write(read)
        except (KeyError):
            continue
    #infile.close()
    #outfile.close()

def index(outfile,edge, I, AF):
    print("index")
    index = "output/bam/coloredBAM_%s_%s_%s.bai" % (edge, I, AF)
    pysam.samtools.index(outfile, index)

infile = pysam.AlignmentFile(bam, "rb")
outfile = pysam.AlignmentFile("output/bam/coloredBAM_NA.bam", "wb", template=infile)

for edge in edges:
    print(edge)
    try:
        write_bam(bam, edge, I, AF)
        print("done")
    except (FileNotFoundError):
        continue


infile.close()
outfile.close()

pysam.samtools.index("output/bam/coloredBAM_NA.bam","output/bam/coloredBAM_NA.bai")