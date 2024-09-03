from os.path import dirname
from sys import path

path.insert( 0 , dirname( __file__ ) ) 

from connector import mass_connector

# https://stackoverflow.com/questions/9427037/relative-path-not-working-even-with-init-py