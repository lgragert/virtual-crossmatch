from rest_framework import serializers


class UNOSSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_UNOS_Antigen_equivalents = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)


	
class AlleleSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_HLA_Alleles = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)
