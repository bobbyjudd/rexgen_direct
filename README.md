# rexgen_direct
Template-free prediction of organic reaction outcomes using graph convolutional neural networks

Described in [A graph-convolutional neural network model for the prediction of chemical reactivity](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04228d)

# Dependencies
- Python (trained/tested using 2.7.6, visualization/deployment compatible with 3.6.1)
- Numpy (trained/tested using 1.12.0, visualization/deployment compatible with 1.14.0)
- Tensorflow (trained/tested using 1.3.0, visualization/deployment compatible with 1.6.0)
- RDKit (trained/tested using 2017.09.1, visualization/deployment compatible with 2017.09.3)
- Django (visualization compatible with 2.0.6)

_note: there may be some issues with relative imports when using Python 2 now; this should be easy to resolve by removing the periods preceding package names_

# Instructions 


### Looking at predictions from the test set
```cd``` into the ```website``` folder and start the Django app using ```python manage.py runserver```. Go to ```http://localhost:8000/visualize``` in a browser to use the interactive visualization tool

### Using the trained models
You can use the fully trained model to predict outcomes by following the example at the end of ```rexgen_direct/rank_diff_wln/directcandranker.py```

### Retraining the models
Look at the two text files in ```rexgen_direct/core_wln_global/notes.txt``` and ```rexgen_direct/rank_diff_wln/notes.txt``` for the exact commands used for training, validation, and testing. You will have to unarchive the data files after cloning this repo.

# Python 3 Instructions
Copy ```1976_Sep2016_USPTOgrants_smiles.rsmi``` into ```rexgen_direct/data```.

Run data preprocessing script

    cd rexgen_direct/data
    python prep_data.py

Create cbond_detailed file

    cd ../core_wln_global

    python nntest_direct.py --test ../data/custom_filtered.rsmi.proc --hidden 300 --depth 3 --model model-300-3-direct --checkpoint ckpt-140000 --verbose 1 --detailed 1 > model-300-3-direct/new_data.cbond_detailed

Get bond predictions - includes reactivity scores in output

    cd ../rank_diff_wln

    python nntest_direct_useScores.py --test ../data/custom_filtered.rsmi.proc --cand ../core_wln_global/model-300-3-direct/new_data.cbond_detailed --hidden 500 --depth 3 --ncand 1500   --ncore 16 --model model-core16-500-3-max150-direct-useScores --checkpoint ckpt-2400000 --verbose 1 > model-core16-500-3-max150-direct-useScores/new_data.cbond_detailed_2400000
    
    python ../scripts/eval_by_smiles.py --gold ../data/1976_Sep2016_USPTOgrants_smiles.rsmi.proc --pred model-core16-500-3-max150-direct-useScores/new_data.cbond_detailed_2400000 --bonds_as_doubles true

Go to website directory and start the Django app:

    cd ../website
    python manage.py runserver


Go to ```http://localhost:8000/visualize/custom``` in a browser to use the interactive visualization tool on the custom dataset predictions
