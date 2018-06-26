
from django.shortcuts import render
from rest_framework.schemas import SchemaGenerator
from rest_framework.permissions import AllowAny

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status, generics
import vxm_hla
import conversion_functions_for_VXM
#from conversion_functions import convert_allele_to_ag
from vxm_hla import allele_truncate
from . import serializers
#from rest_framework_swagger.views import get_swagger_view

import virtual_crossmatch

from virtual_crossmatch import vxm_allele_codes


pop_acro_dict = {"AFA": "African American", "API": "Asia/Pacific Islander", "CAU": "Caucasian", "HIS": "Hispanic", 
"NAM": "Native American", "AAFA": "African American", "AFB": "African Black",  "AINDI": "South Asian Indian", 
"AISC": "American Indian-South or Central American", "ALANAM": "Alaska Native", "AMIND": "North American Indian", 
"CARB": "Caribbean Black", "CARHIS": "Caribbean Hispanic", "CARIBI": "Caribbean Indian", "EURACU": "European Caucasian", 
"FILII": "Filipino", "HAWI": "Hawaiian or Pacific Islander", "JAPI": "Japanese", "KORI": "Korean", 
"MENAFC": "Middle Eastern or N. Coast of Africa", "MSWHIS": "Mexican or Chicano HIS", "NCHI": "Chinese", 
"SCAHIS": "South or Central AMerican Hispanic", "SCAMB": "South or Central American Black", 
"SCSEAI": "South East Asian", "VIET": "Vietnamese"}



class MultipleAlleleCodesApiView(generics.GenericAPIView):
    serializer_class = serializers.VICTOR_MACSerializer

    def post(self, request, format=None):
        """Returns UNOS antigen for an allele."""
        """parameters: 
        allele:string
        """
        serializer = serializers.VICTOR_MACSerializer(data=request.data)    

        if serializer.is_valid(raise_exception=True):
            
            #serializer.create(allele)
            #allele = serializer.save()
            #serializer.create()
            donor_macs = serializer.data.get('Donor_Allele_Codes')
            donor_macs_list_split = donor_macs.split(" ")
            pop = serializer.data.get('Donor_ethinicity')
            ua_list = serializer.data.get('Candidate_Unacceptable_antigens')
            ua_list_split = ua_list.split(" ")
            output = virtual_crossmatch.vxm_allele_codes(donor_macs_list_split, pop, ua_list_split)
            popSpecFul = pop_acro_dict[pop]
            Donor_ags = ", ".join(output[0])
            Candidate_unacceptable_ags = ", ".join(output[1])
            print(output)
            if len(output[2]) != 0:

                VXM_out = "Positive"
                conflicts =  ", ".join(output[2])
                probs = output[3]

                return Response({"VICTOR VXM": VXM_out, "Donor UNOS Antigen Equivalents": Donor_ags, "Donor Race": popSpecFul,"Candidate Unacceptable Antigens": Candidate_unacceptable_ags, 
                    "Conflicting Antigens": conflicts, "Probabilities of virtual crossmatch": probs },  status=status.HTTP_200_OK)

            else: 
                VXM_out = "Negative"
                return Response({"VICTOR VXM": VXM_out, "Donor UNOS Antigen Equivalents": Donor_ags, "Donor Race": popSpecFul, "Candidate Unacceptable Antigens": Candidate_unacceptable_ags}, 
                 status=status.HTTP_200_OK)


        else:
            return Response({"Error": "Check if allele is IMGT/HLA"}, status=status.HTTP_400_BAD_REQUEST)
            #serializer.errors, status=status.HTTP_400_BAD_REQUEST)    