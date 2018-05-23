
from django.contrib import admin
from django.urls import path
from vxm_app import views_vxm

urlpatterns = [
    path('', views_vxm.vxm_home, name="victor_home"),
    path('glsvxm/', views_vxm.gl_string, name="glsvxm"),
    path('matchgls/', views_vxm.match_gl, name="matchgls"),
    path('macvxm/', views_vxm.multiple_allele_codes, name="macvxm"),
    path('matchmac/', views_vxm.match_ac, name="matchmac"),
    path('admin/', admin.site.urls),
    path('VICTORLICENSE/', views_vxm.victor_license, name="victor_license"),

]
