from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, Http404
from .data import do_index
import base64
from io import BytesIO
import json 
from .data import df_rank_validation, df_rank_test

def index(request, dataset='default'):
    print('index: {}'.format(dataset))

    if dataset not in set(['default', 'test', 'validation']):
        raise Http404("No dataset matches that query.")

    return render(request, 'example.html', {'dataset':dataset})

def display(request):
    print('Got request')

    data = {'err': False}
    dataset = request.GET.get('dataset')
    print('display: {}'.format(dataset))
    index = int(request.GET.get('index', 1))
    showmap = json.loads(request.GET.get('showmap', 'true'))
    highlight = int(request.GET.get('highlight', 0))

    print(index)
    atts = [highlight] if highlight else []
    if dataset == 'validation':
        index = min(index, len(df_rank_validation))
        entry,img = do_index(index - 1, showmap=showmap, atts=atts, df_rankpred=df_rank_validation)
    elif dataset == 'test':
        index = min(index, len(df_rank_test))
        entry,img = do_index(index - 1, showmap=showmap, atts=atts, df_rankpred=df_rank_test)
    else:
        entry,img = do_index(index - 1, showmap=showmap, atts=atts)
    buffered = BytesIO()
    img.save(buffered, format="png")
    img_str = base64.b64encode(buffered.getvalue()).decode()

    data['img_str'] = img_str
    data['index'] = index
    data['reactants'] = list(entry['reactants'])
    data['products'] = list(entry['products'])

    return JsonResponse(data)