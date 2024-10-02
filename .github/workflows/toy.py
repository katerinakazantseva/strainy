import gfapy

g2 = gfapy.Gfa.from_file("out_strainy/strain_contigs.gfa")

def check_gfa(path):
    g2 = gfapy.Gfa.from_file(path)
    return [len(g2.segment_names),len(g2.edges)]


def test_answer():
    assert check_gfa("out_strainy/strain_contigs.gfa") ==[15,16]
    assert check_gfa("out_strainy/strain_unitigs.gfa") ==[22,23]
