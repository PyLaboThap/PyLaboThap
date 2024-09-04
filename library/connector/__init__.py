from os.path import dirname
from sys import path

path.insert( 0 , dirname( __file__ ) ) 

# from connector import mass_connector

# https://stackoverflow.com/questions/9427037/relative-path-not-working-even-with-init-py

# library/connector/__init__.py

# Only include necessary imports
from .mass_connector import MassConnector
from .work_connector import WorkConnector
from .heat_connector import HeatConnector