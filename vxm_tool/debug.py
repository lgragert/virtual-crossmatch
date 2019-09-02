#! usr/bin/python

from reverse_conversion import map_ag_for_proposed_algo

from vxm_hla import ags_to_strings

from conversion_functions_for_VXM import gl_string_ags, convert_ag_list_to_gls

ag1 = "DQ6"
ag2 = "DQ2"

alleles1 = map_ag_for_proposed_algo(ag1)
alleles2 = map_ag_for_proposed_algo(ag2)

string = ags_to_strings(ag1, ag2)
#print(len(string))

# no of alleles of DQ6 : 258
# no of alleles of DQ2 : 101
# no of possible genotypes = 26058
pop = "CAU"
x = gl_string_ags(string, pop)
print(x[1][0])