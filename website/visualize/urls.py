from django.urls import path

from . import views

urlpatterns = [
    path('ajax', views.display, name='display'),
    path('', views.index, name='index'),
    path('<dataset>', views.index, name='index')
    
]