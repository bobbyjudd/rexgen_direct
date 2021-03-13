from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from .data import do_index
import base64
from io import BytesIO
import json 
from .data import df_rank_custom

def index(request):
    return render(request, 'example.html', {})

def index_custom(request):
    return render(request, 'custom_example.html', {})

def display(request):
    print('Got request')

    data = {'err': False}
    index = int(request.GET.get('index', 1))
    showmap = json.loads(request.GET.get('showmap', 'true'))
    highlight = int(request.GET.get('highlight', 0))

    print(index)
    atts = [highlight] if highlight else []
    img = do_index(index - 1, showmap=showmap, atts=atts)
    buffered = BytesIO()
    img.save(buffered, format="png")
    img_str = base64.b64encode(buffered.getvalue()).decode()

    data['img_str'] = img_str
    data['index'] = index

    return JsonResponse(data)

def display_custom(request):
    print('Got request')

    data = {'err': False}
    index = int(request.GET.get('index', 1))
    showmap = json.loads(request.GET.get('showmap', 'true'))
    highlight = int(request.GET.get('highlight', 0))

    print(index)
    atts = [highlight] if highlight else []
    img = do_index(index - 1, showmap=showmap, atts=atts, df_rankpred=df_rank_custom)
    buffered = BytesIO()
    img.save(buffered, format="png")
    img_str = base64.b64encode(buffered.getvalue()).decode()

    data['img_str'] = img_str
    data['index'] = index

    return JsonResponse(data)