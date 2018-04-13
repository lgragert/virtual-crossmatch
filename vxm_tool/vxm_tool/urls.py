"""vxm_tool URL Configuration

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
from django.urls import path
from vxm_app import views

urlpatterns = [
    path('', views.vxm_home, name="home"),
    path('glsvxm/', views.gl_string, name="glsvxm"),
    path('matchgls/', views.match_gl, name="matchgls"),
    path('macvxm/', views.multiple_allele_codes, name="macvxm"),
    path('matchmac/', views.match_ac, name="matchmac"),
    path('admin/', admin.site.urls),
    path('license/', views.license, name="license"),

]
