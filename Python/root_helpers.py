def print_daughters(mdst, i_part, indent=0):
    """Given an index to an MC particle, print all the daughters.
    This function will recursively print the daughters of the daughters.
    Arguments:
        mdst   -- a MiniDst object that was linked to a TTree
        i_part -- the index of the particle to print
        ident  -- the amount of indentation of the output.
    """

    print(" "*indent+f" {i_part:3d}  pdg: {mdst.mc_part_pdg_id[i_part]:4d}  E: {mdst.mc_part_energy[i_part]:9.6f} " +
    f"p = ({mdst.mc_part_px[i_part]:9.6f},{mdst.mc_part_py[i_part]:9.6f},{mdst.mc_part_pz[i_part]:9.6f})"          )
    if len(mdst.mc_part_daughters[i_part]) > 0:
        print(" "*(indent+14) + "| ")
        for i in range(len(mdst.mc_part_daughters[i_part])):
            ii = mdst.mc_part_daughters[i_part][i]  # Get the daughter reference
            print_daughters(mdst, ii, indent+11)            # Print by recursing

def print_mc_particle_tree(mdst):
    """Print the MCParticle tree.
    Arguments:
        mdst -- a MiniDst object that was linked to a TTree.
    """
    for i in range(len(mdst.mc_part_parents)):
        if len(mdst.mc_part_parents[i]) == 0:  # top level particle
            print_daughters(mdst, i, 0)


def get_vertex_dictionary():
    """Return a dictionary that translates the Vertex numbers to the
    ParticleType Enum name"""

    out = {
        0: "FINAL_STATE_PARTICLE_KF",
        1: "UC_V0_VERTICES_KF",
        2: "BSC_V0_VERTICES_KF",
        3: "TC_V0_VERTICES_KF",
        4: "UC_MOLLER_VERTICES_KF",
        5: "BSC_MOLLER_VERTICES_KF",
        6: "TC_MOLLER_VERTICES_KF",
        7: "OTHER_ELECTRONS_KF",
        8: "UC_VC_VERTICES_KF",
        0 + 9: "FINAL_STATE_PARTICLE_GBL",
        1 + 9: "UC_V0_VERTICES_GBL",
        2 + 9: "BSC_V0_VERTICES_GBL",
        3 + 9: "TC_V0_VERTICES_GBL",
        4 + 9: "UC_MOLLER_VERTICES_GBL",
        5 + 9: "BSC_MOLLER_VERTICES_GBL",
        6 + 9: "TC_MOLLER_VERTICES_GBL",
        7 + 9: "OTHER_ELECTRONS_GBL",
        8 + 9: "UC_VC_VERTICES_GBL",
    }
    return out


def SetStyle(choice=0):
    if "R" not in globals():
        import ROOT as R
    if "gROOT" not in globals():
        from ROOT import gROOT

    hpsStyle = R.TStyle("HPS", "HPS style")

    # use plain black on white colors
    icol = 0
    hpsStyle.SetFrameBorderMode(icol)
    hpsStyle.SetCanvasBorderMode(icol)
    hpsStyle.SetPadBorderMode(icol)
    hpsStyle.SetPadColor(icol)
    hpsStyle.SetCanvasColor(icol)
    hpsStyle.SetStatColor(icol)
    # hpsStyle.SetFillColor(icol)

    # set the paper & margin sizes
    hpsStyle.SetPaperSize(20,26)
    hpsStyle.SetPadTopMargin(0.05)
    hpsStyle.SetPadRightMargin(0.05)
    hpsStyle.SetPadBottomMargin(0.18)
    hpsStyle.SetPadLeftMargin(0.14)

    # use large fonts
    # font=72
    if(choice == 0):
        font = 42         # helvetica-medium-r-normal
        title_font = 42
        title_size = 0.08
        title_size_z = 0.045

        label_font = 42
        label_size = 0.08
        label_size_z = 0.045
        hpsStyle.SetOptTitle(0)
        hpsStyle.SetOptStat(0)
        hpsStyle.SetOptFit(0)

    elif(choice == 1):
        #font=72
        font = 132        # times-medium-r-normal
        title_font = 132
        title_size = 0.08
        title_size_z = 0.045

        label_font = 42
        label_size = 0.035
        label_size_z = 0.035
        hpsStyle.SetOptTitle(0)
        hpsStyle.SetOptStat(0)
        hpsStyle.SetOptFit(0)


    hpsStyle.SetTextFont(font)
    hpsStyle.SetTextSize(label_size)

    hpsStyle.SetLabelFont(label_font, "x")
    hpsStyle.SetTitleFont(title_font, "x")
    hpsStyle.SetLabelFont(label_font, "y")
    hpsStyle.SetTitleFont(title_font, "y")
    hpsStyle.SetLabelFont(label_font, "z")
    hpsStyle.SetTitleFont(title_font, "z")

    hpsStyle.SetLabelSize(label_size, "x")
    hpsStyle.SetTitleSize(title_size, "x")
    hpsStyle.SetLabelSize(label_size, "y")
    hpsStyle.SetTitleSize(title_size, "y")
    hpsStyle.SetLabelSize(label_size_z, "z")
    hpsStyle.SetTitleSize(title_size_z, "z")

    hpsStyle.SetTitleOffset(0.7, "y")
    hpsStyle.SetTitleOffset(1.15, "x")

    #use bold lines and markers
    #hpsStyle.SetMarkerStyle(20)
    hpsStyle.SetMarkerSize(1.0)
    hpsStyle.SetHistLineWidth(1)
    hpsStyle.SetLineStyleString(2, "[12 12]")  # postscript dashes

    #get rid of X error bars and y error bar caps
    #hpsStyle.SetErrorX(0.001)

    #do not display any of the standard histogram decorations

    # put tick marks on top and RHS of plots
    hpsStyle.SetPadTickX(1)
    hpsStyle.SetPadTickY(1)

    #gROOT.SetStyle("Plain")
    gROOT.SetStyle("HPS")
    gROOT.ForceStyle()

    NCont = 255
    stops = R.std.vector("double")([0.00, 0.34, 0.61, 0.84, 1.00])
    red = R.std.vector("double")([0.00, 0.00, 0.87, 1.00, 0.51])
    green = R.std.vector("double")([0.00, 0.81, 1.00, 0.20, 0.00])
    blue = R.std.vector("double")([0.51, 1.00, 0.12, 0.00, 0.00])
    R.TColor.CreateGradientColorTable(stops.size(), stops.data(), red.data(), green.data(), blue.data(), NCont)
    R.gStyle.SetNumberContours(NCont)
