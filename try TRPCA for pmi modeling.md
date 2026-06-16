
https://github.com/tydymy/TRPCA
```
conda create --name torch_env

conda activate torch_env
conda install pytorch torchvision torchaudio -c pytorch
#verify
python
import torch  
  
print(torch.__version__)  
print(torch.cuda.is_available())
```


install dependencies within that env

```
module load qiime2/2026.1_amplicon # Core dependencies
pip install torch numpy pandas scikit-learn

# Microbiome analysis
pip install biom-format
pip install cython
pip install "scikit-bio==0.5.9"  
pip install gemelli

# Visualization and analysis
pip install matplotlib seaborn shap optuna tqdm
```

## Installation

1. Clone the repository:

```shell
git clone https://github.com/tydymy/TRPCA.git
cd TRPCA
git checkout trpca_v1
```

```
pip install jupyter
```

```
pip install "scikit-bio==0.5.9"
```
this is becoming too crazy to figure out all the matching dependencies etc. taking a break from this for right now