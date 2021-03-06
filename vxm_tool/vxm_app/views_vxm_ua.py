
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

from virtual_crossmatch import vxm_uags

class UNOSAgsApiView(generics.GenericAPIView):
    serializer_class = serializers.VICTOR_UNOSSerializer

    def post(self, request, format=None):
        """Computes VXM for UNOS Ags"""
        """parameters: 
        allele:string
        """
        serializer = serializers.VICTOR_UNOSSerializer(data=request.data)    

        if serializer.is_valid(raise_exception=True):
            
            #serializer.create(allele)
            #allele = serializer.save()
            #serializer.create()
            unos_ags_list = serializer.data.get('Donor_UNOS_Antigen_equivalents')
            #print(unos_ags_list)
            ua_list = serializer.data.get('Candidate_Unacceptable_antigens')
            unos_ag_split = unos_ags_list.split(" ")
            #print(unos_ag_split)
            ua_list_split = ua_list.split(" ")
            #print(ua_list_split)
            
            output = virtual_crossmatch.vxm_uags(unos_ag_split, ua_list_split)
            Donor_ags = ", ".join(output[0])
            Candidate_unacceptable_ags = ", ".join(output[1])
            #print(len(output))
            if len(output[2]) != 0:

                VXM_out = "Positive"
                conflicts =  ", ".join(output[2])

                return Response({"VICTOR VXM": VXM_out, "Donor UNOS Antigen Equivalents": Donor_ags, "Candidate Unacceptable Antigens": Candidate_unacceptable_ags, "Conflicting Antigens": conflicts, },  status=status.HTTP_200_OK)

            else: 
                VXM_out = "Negative"
                return Response({"VICTOR VXM": VXM_out, "Donor UNOS Antigen Equivalents": Donor_ags, "Candidate Unacceptable Antigens": Candidate_unacceptable_ags}, 
                 status=status.HTTP_200_OK)


        else:
            return Response({"Error": "Check if allele is IMGT/HLA"}, status=status.HTTP_400_BAD_REQUEST)
            #serializer.errors, status=status.HTTP_400_BAD_REQUEST)    