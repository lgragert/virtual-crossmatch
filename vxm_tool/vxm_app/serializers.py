from rest_framework import serializers


class UNOSSerializer(serializers.Serializer):
	"""Serializes a name field for testing our APIView."""
	unos_ags_list = serializers.CharField(max_length=None)
	ua_list = serializers.CharField(max_length=None)
