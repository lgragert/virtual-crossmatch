from rest_framework import serializers


class VICTOR_UNOSSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_UNOS_Antigen_equivalents = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)


class VICTOR_Prop_UNOSSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_UNOS_Antigen_equivalents = serializers.CharField(max_length=None)
	Donor_ethinicity = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)


	
class VICTOR_AlleleSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_HLA_Alleles = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)


class VICTOR_GLstringSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_GL_string = serializers.CharField(max_length=None)
	Donor_ethinicity = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)
	

class VICTOR_MACSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_Allele_Codes = serializers.CharField(max_length=None)
	Donor_ethinicity = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)
