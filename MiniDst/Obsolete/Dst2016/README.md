# Dst2016

The Dst2016 class derives from the BaseAna class and exists only to allow you
to read the old 2016 style DST files and convert them to the flat TTree ntuples that
are used by the rest of the code here. The executable make_mini_dst will use this
class when presented with ROOT file input to convert it to the flat ROOT file output.
