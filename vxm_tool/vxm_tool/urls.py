
from django.contrib import admin
from django.urls import path
from django.conf import settings
from vxm_app import views_vxm, views_vxm_ua, views_vxm_alleles, views_vxm_gls

urlpatterns = [
    path('', views_vxm.vxm_home, name="victor_home"),
    path('unosags/', views_vxm.unos_ags, name="unosags"),
    path('matchuags/', views_vxm.match_ags, name="uagsvxm"),
    path('highresallele/', views_vxm.highres_allele, name="highresallele"),
    path('matchhiresallele/', views_vxm.match_hi_res_alleles, name="matchallele"),
    path('glsvxm/', views_vxm.gl_string, name="glsvxm"),
    path('matchgls/', views_vxm.match_gl, name="matchgls"),
    path('macvxm/', views_vxm.multiple_allele_codes, name="macvxm"),
    path('matchmac/', views_vxm.match_ac, name="matchmac"),
    path('admin/', admin.site.urls),
    path('VICTORLICENSE/', views_vxm.victor_license, name="victor_license"),
    path('ags_unos/', views_vxm_ua.UNOSAgsApiView.as_view()),
    path('victor_alleles/', views_vxm_alleles.AllelesApiView.as_view()),
    path('victor_gls/', views_vxm_gls.GenotypeListStringApiView.as_view()),



]
