"""transplanttoolbox URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf.urls import url
from django.conf import settings
from vxm_app import views_vxm, views_vxm_ua, views_vxm_alleles, views_vxm_gls, views_vxm_macs
from rest_framework_swagger.views import get_swagger_view

schema_view = get_swagger_view(title='Web services interface information for Transplanttoolbox')

from django.urls import path

urlpatterns = [
    path('', views_vxm.vxm_home, name="victor_home"),
    path('services', schema_view, name="victor_services"),
    path('unosags/', views_vxm.unos_ags, name="unosags"),
    path('matchuags/', views_vxm.match_ags, name="uagsvxm"),
    path('highresallele/', views_vxm.highres_allele, name="highresallele"),
    path('matchhiresallele/', views_vxm.match_hi_res_alleles, name="matchallele"),
    path('proposeduags/', views_vxm.proposed_unos_ags, name="proposeduags"),
    path('proposedmatchuags/', views_vxm.match_proposed_uags, name="proposedmatchuags"),
    path('glsvxm/', views_vxm.gl_string, name="glsvxm"),
    path('matchgls/', views_vxm.match_gl, name="matchgls"),
    path('glsextendedtable/', views_vxm.gl_string_extended_table, name="gl_string_extended_table"),
    path('macvxm/', views_vxm.multiple_allele_codes, name="macvxm"),
    path('matchmac/', views_vxm.match_ac, name="matchmac"),
    path('admin/', admin.site.urls),
    path('VICTORLICENSE/', views_vxm.victor_license, name="victor_license"),
    path('ags_unos/', views_vxm_ua.UNOSAgsApiView.as_view()),
    path('highres_alleles/', views_vxm_alleles.AllelesApiView.as_view()),
    path('gls/', views_vxm_gls.GenotypeListStringApiView.as_view()),
    path('macs/', views_vxm_macs.MultipleAlleleCodesApiView.as_view()),


]
