from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('ajax', views.display, name='display'),
    path('custom', views.index_custom, name='index_custom'),
    path('customajax', views.display_custom, name='display_custom')
]