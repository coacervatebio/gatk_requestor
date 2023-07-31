import yaml
from pathlib import Path

ROOT_PATH = Path(__file__).resolve().parent
CONFIG_PATH = ROOT_PATH.joinpath('config.yaml')

with open(CONFIG_PATH, 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

print(config)