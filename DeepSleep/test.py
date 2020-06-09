import sys
import pandas as pd 
# testing pbs batch system
print(pd.__version__)
print(f"Hello World {sys.argv[1]}")
x = pd.DataFrame({'col1':[1,2,3]})
x.to_pickle(f'{sys.argv[1]}_test.pkl')
