from rest_framework import serializers


class UNOSSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	Donor_UNOS_Antigen_equivalents = serializers.CharField(max_length=None)
	Candidate_Unacceptable_antigens = serializers.CharField(max_length=None)
